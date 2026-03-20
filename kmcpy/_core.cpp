/*
 * kmcpy/_core.cpp
 *
 * pybind11 binding between Python and libkmc_core.a.
 *
 * Stage 1: bins are kept in RAM while max_ram_gb allows.
 *          When exceeded, KMC automatically spills bins to tmp_path.
 *          ramOnlyMode=false — KMC decides bin-by-bin whether to spill.
 *
 * Stage 2: output (packedKmers/prefixArray) always comes back to RAM.
 *          withoutOutput=true — no .kmc_pre/.kmc_suf written to disk.
 *          strict_mem=true (default) — Stage 2 sorting never exceeds RAM.
 *
 * ============================================================================
 * OUTPUT ARRAYS  (all built in one pass)
 * ============================================================================
 *
 * kmers  (n, k)  uint8
 *   One cell per base position.  Value = raw KMC 2-bit code:
 *     A=0  C=1  G=2  T=3
 *   Zero transformation from the bin format.
 *   Valid for ALL k (1–256).
 *   Used by decoded=True path in Python to produce ACGT strings.
 *   Always allocated regardless of drop_count.
 *
 * key_words  (n, n_words)  uint64
 *   Space-optimal packed k-mer key for decoded=False.
 *   n_words = ceil(k * 2 / 64) — the minimum number of uint64 words
 *   required to hold 2*k bits with no wasted words (at most 63 wasted bits).
 *
 *   Bit layout within the n_words * 64-bit block:
 *     word[0] holds the MOST SIGNIFICANT bits (base[0] side).
 *     word[n_words-1] holds the LEAST SIGNIFICANT bits (base[k-1] side).
 *     The key is RIGHT-ALIGNED within the word block:
 *       base[0]   occupies bits [2k-1 : 2k-2] counting from bit 0 of word[0]
 *       base[k-1] occupies bits [1    : 0    ] of word[n_words-1]
 *       Unused HIGH bits of word[0] are always zero.
 *     A=00, C=01, G=10, T=11  (exact KMC 2-bit encoding, preserved as-is).
 *
 *   Space efficiency by k range:
 *     k ≤  32: n_words=1  — 1 × uint64 column  "kmer_0"
 *     k ≤  64: n_words=2  — 2 × uint64 columns "kmer_0", "kmer_1"
 *     k ≤  96: n_words=3  — 3 × uint64 columns
 *     k ≤ 128: n_words=4  — 4 × uint64 columns
 *     k ≤ 192: n_words=6  — 6 × uint64 columns
 *     k ≤ 256: n_words=8  — 8 × uint64 columns (KMC maximum)
 *
 *   cuDF / pandas merge key:
 *     cols = [f"kmer_{i}" for i in range(n_words)]
 *     merged = gdf1.merge(gdf2, on=cols)
 *
 * counts  (n,)  uint32   — occurrence counts, decoded from little-endian bin bytes.
 *   ONLY allocated and transferred when drop_count=false (default).
 *   When drop_count=true, this array is never allocated — saving n*4 bytes
 *   of heap allocation, fill time, and Python/C++ transfer overhead.
 *   For 100M k-mers: 400 MB saved.  For 8 parallel jobs: 3.2 GB saved.
 *
 * ============================================================================
 * IMPLEMENTATION NOTES
 * ============================================================================
 *
 * KMC invariants we rely on:
 *   1. pfx_len ≤ 15  →  prefix always fits in a uint64 (≤ 30 bits used).
 *   2. suffix_len % 4 == 0  →  suf_bytes = suf_bases / 4 exactly, no padding.
 *   3. Each suffix byte holds exactly 4 bases, packed MSB-first.
 *
 * Key assembly uses a stack-allocated uint64[8] buffer (512 bits).
 * No heap allocation per k-mer.  No __uint128_t dependency.
 * The buffer is always zeroed and filled from MSB to LSB to match the
 * right-aligned layout described above.
 *
 * When drop_count=true, count bytes in packed_buf are skipped entirely —
 * we still advance rec by stride (suf_bytes + cnt_bytes) to reach the next
 * record, but we never read the count bytes or write to counts_arr.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <cstring>   // memset
#include <cstdint>

// Fork: kmc_runner.h is copied to CMAKE_BINARY_DIR/include/ at configure time.
// ${CMAKE_BINARY_DIR}/include is on the include path -- so just "kmc_runner.h".
#include "kmc_runner.h"

namespace py = pybind11;


// ---------------------------------------------------------------------------
// KmerKeyBuffer
//
// Stack-allocated 512-bit (8 × uint64) buffer for assembling the packed
// k-mer key.  Indexed word[0] = most significant, word[7] = least significant.
//
// Only the top n_words words are populated; the rest are always zero.
// After assembly, word[0..n_words-1] are copied into the numpy output array.
// ---------------------------------------------------------------------------
struct KmerKeyBuffer
{
    uint64_t word[8];   // 512 bits — covers k up to 256

    void clear() { std::memset(word, 0, sizeof(word)); }

    // Set a single 2-bit code at base position p (0 = leftmost base).
    // total_bits = k * 2.
    // The key is right-aligned: base[k-1] is at bit 1:0 of word[n_words-1].
    //
    // Simplified with right-alignment (total field = n_words * 64):
    //   field_bit = total_bits - 2 - 2*p   (LSB of the 2-bit code)
    //   word_idx  = n_words - 1 - (field_bit / 64)
    //   shift     = field_bit % 64
    inline void set_base(uint32_t p, uint8_t code,
                         uint32_t total_bits, uint32_t n_words) noexcept
    {
        uint32_t field_bit = total_bits - 2u - 2u * p;
        uint32_t widx      = n_words - 1u - (field_bit / 64u);
        uint32_t shift     = field_bit % 64u;
        word[widx] |= static_cast<uint64_t>(code) << shift;
    }
};


// ---------------------------------------------------------------------------
// _fill_key_buffer
//
// Populate a KmerKeyBuffer from KMC's split prefix + suffix representation.
//
// Parameters:
//   buf        — output buffer (must be cleared before call)
//   pfx        — uint64 LUT row index; encodes pfx_len bases, 2 bits each,
//                MSB = base[0].  Max pfx_len=15 → max 30 bits → fits in uint64.
//   rec        — pointer into packed_buf; first suf_bytes bytes are the suffix.
//                KMC packs 4 bases per byte, MSB-first.
//                suffix_len % 4 == 0 always (KMC invariant) → no partial bytes.
//   pfx_len    — number of prefix bases (0–15).
//   suf_bytes  — number of suffix bytes = suf_bases / 4.
//   suf_bases  — kmer_len - pfx_len.
//   total_bits — kmer_len * 2.
//   n_words    — ceil(total_bits / 64).
// ---------------------------------------------------------------------------
static void _fill_key_buffer(
    KmerKeyBuffer&  buf,
    uint64_t        pfx,
    const uint8_t*  rec,
    uint32_t        pfx_len,
    uint32_t        suf_bytes,
    uint32_t        suf_bases,
    uint32_t        total_bits,
    uint32_t        n_words) noexcept
{
    // --- Prefix bases (0 .. pfx_len-1) ---
    // pfx holds pfx_len*2 bits in its LSBs, MSB = base[0].
    for (uint32_t p = 0; p < pfx_len; ++p)
    {
        uint32_t shift = 2u * (pfx_len - 1u - p);
        uint8_t  code  = static_cast<uint8_t>((pfx >> shift) & 0x3u);
        buf.set_base(p, code, total_bits, n_words);
    }

    // --- Suffix bases (pfx_len .. kmer_len-1) ---
    // Each byte holds 4 bases, MSB-first.
    // suffix_len % 4 == 0 → each byte is always full (4 bases).
    uint32_t base_idx = pfx_len;
    for (uint32_t b = 0; b < suf_bytes; ++b)
    {
        uint8_t bval = rec[b];
        buf.set_base(base_idx,     (bval >> 6) & 0x3u, total_bits, n_words);
        buf.set_base(base_idx + 1, (bval >> 4) & 0x3u, total_bits, n_words);
        buf.set_base(base_idx + 2, (bval >> 2) & 0x3u, total_bits, n_words);
        buf.set_base(base_idx + 3,  bval        & 0x3u, total_bits, n_words);
        base_idx += 4u;
    }
}


static py::dict _count_kmers_internal(
    const std::vector<std::string>& input_files,
    const std::string&  tmp_path,
    uint32_t    kmer_len,
    uint32_t    max_ram_gb,
    uint32_t    n_threads,
    uint32_t    n_readers,
    uint32_t    n_splitters,
    uint32_t    signature_len,
    uint32_t    n_bins,
    bool        canonical,
    bool        homopolymer,
    bool        strict_mem,
    int         input_type,
    uint64_t    cutoff_min,
    uint64_t    cutoff_max,
    uint32_t    counter_max,
    bool        drop_count,      // Fork: if true, skip counts_arr entirely
    bool        need_bases       // Fork: if false, skip kmers_arr (decoded=False path)
)
{
    KMC::Stage1Params s1;
    s1.SetInputFiles(input_files)
      .SetKmerLen(kmer_len)
      .SetRamOnlyMode(false)
      .SetTmpPath(tmp_path)
      .SetCanonicalKmers(canonical)
      .SetHomopolymerCompressed(homopolymer)
      .SetMaxRamGB(max_ram_gb)
      .SetNThreads(n_threads)
      .SetSignatureLen(signature_len);

    if (n_readers   > 0) s1.SetNReaders(n_readers);
    if (n_splitters > 0) s1.SetNSplitters(n_splitters);
    if (n_bins      > 0) s1.SetNBins(n_bins);

    switch (input_type) {
        case 1:  s1.SetInputFileType(KMC::InputFileType::FASTA);           break;
        case 2:  s1.SetInputFileType(KMC::InputFileType::MULTILINE_FASTA); break;
        case 3:  s1.SetInputFileType(KMC::InputFileType::BAM);             break;
        case 4:  s1.SetInputFileType(KMC::InputFileType::KMC);             break;
        default: s1.SetInputFileType(KMC::InputFileType::FASTQ);           break;
    }

    KMC::Stage2Params s2;
    s2.SetCutoffMin(cutoff_min)
      .SetCutoffMax(cutoff_max)
      .SetCounterMax(counter_max)
      .SetMaxRamGB(max_ram_gb)
      .SetNThreads(n_threads)
      .SetWithoutOutput(true)
      .SetStrictMemoryMode(strict_mem);

    KMC::Runner runner;
    KMC::Stage1Results r1;
    KMC::Stage2Results r2;

    {
        py::gil_scoped_release _release;
        r1 = runner.RunStage1(s1);
        r2 = runner.RunStage2(s2);
    }

    const auto&    packed     = r2.packedKmers;
    const auto&    prefixes   = r2.prefixArray;
    const size_t   n          = prefixes.size();
    const uint32_t suf_bytes  = r2.kmerSufBytes;
    const uint32_t cnt_bytes  = r2.counterSize;
    const uint32_t pfx_len    = r2.lutPrefixLen;
    const uint32_t stride     = suf_bytes + cnt_bytes;
    const uint32_t suf_bases  = kmer_len - pfx_len;
    const uint32_t total_bits = kmer_len * 2u;

    // n_words = ceil(k * 2 / 64) — minimum uint64 words to hold all bits.
    const uint32_t n_words = (total_bits + 63u) / 64u;

    // small-k path: k <= 13.
    const bool small_k = (suf_bytes == 0 && pfx_len == 0 && kmer_len <= 13);

    // -----------------------------------------------------------------------
    // Allocate output arrays.
    //
    // kmers_arr  (n, k)        uint8  — per-base 2-bit codes.
    //   ONLY allocated when need_bases=true (decoded=True path).
    //   When need_bases=false, this is a zero-element array — no allocation,
    //   no fill, no transfer. Saves n*k bytes at the C++→Python boundary.
    //   Memory saving: 25 MB per 1M k-mers at k=25.
    //
    // key_words  (n, n_words)  uint64 — packed k-mer key, always built.
    //   Space-optimal: n_words = ceil(k*2/64).
    //   Used directly for decoded=False output AND as source for decoded=True.
    //
    // counts_arr (n,)          uint32 — ONLY allocated when drop_count=false.
    //   Memory saving: 400 MB per 100M k-mers when drop_count=true.
    // -----------------------------------------------------------------------
    py::array_t<uint8_t>  kmers_arr  ({
        need_bases ? (py::ssize_t)n : (py::ssize_t)0,
        need_bases ? (py::ssize_t)kmer_len : (py::ssize_t)0
    });
    py::array_t<uint64_t> key_words  ({ (py::ssize_t)n, (py::ssize_t)n_words  });
    py::array_t<uint32_t> counts_arr ({
        drop_count ? (py::ssize_t)0 : (py::ssize_t)n
    });

    {
        // km accessor — only valid when need_bases=true
        uint8_t* km_ptr = need_bases ? kmers_arr.mutable_data() : nullptr;
        auto kw  = key_words.mutable_unchecked<2>();

        // counts accessor — only used when drop_count=false.
        // Declared outside the loop to avoid repeated branching.
        uint32_t* cnt_ptr = drop_count
            ? nullptr
            : counts_arr.mutable_data();

        KmerKeyBuffer buf;

        for (size_t i = 0; i < n; ++i)
        {
            // ----------------------------------------------------------
            // Small-k path (k <= 13):
            //   prefixes[i] IS the full right-aligned packed integer.
            //   n_words = 1 for all k <= 13 (max 26 bits < 64).
            // ----------------------------------------------------------
            if (small_k)
            {
                kw(i, 0) = prefixes[i];

                // Fill per-base array only when needed (decoded=True)
                if (km_ptr)
                {
                    uint64_t kdata = prefixes[i];
                    uint8_t* row   = km_ptr + i * kmer_len;
                    for (int32_t p = (int32_t)kmer_len - 1; p >= 0; --p)
                    {
                        row[p] = static_cast<uint8_t>(kdata & 0x3u);
                        kdata >>= 2u;
                    }
                }

                // Count — only decode when needed
                if (cnt_ptr)
                {
                    uint32_t count = 0;
                    const uint8_t* rec = packed.data() + i * cnt_bytes;
                    for (uint32_t b = 0; b < cnt_bytes; ++b)
                        count |= static_cast<uint32_t>(rec[b]) << (b * 8u);
                    cnt_ptr[i] = count;
                }
                continue;
            }

            // ----------------------------------------------------------
            // Normal path (k >= 14):
            // ----------------------------------------------------------
            const uint8_t* rec = packed.data() + i * stride;
            uint64_t pfx       = prefixes[i];

            // --- Build packed key (always) ---
            buf.clear();
            _fill_key_buffer(buf, pfx, rec,
                             pfx_len, suf_bytes, suf_bases,
                             total_bits, n_words);
            for (uint32_t w = 0; w < n_words; ++w)
                kw(i, w) = buf.word[w];

            // --- Fill per-base array only when needed (decoded=True) ---
            if (km_ptr)
            {
                uint8_t* row = km_ptr + i * kmer_len;

                for (uint32_t p = 0; p < pfx_len; ++p)
                    row[p] = static_cast<uint8_t>(
                        (pfx >> (2u * (pfx_len - 1u - p))) & 0x3u);

                uint32_t base_idx = pfx_len;
                for (uint32_t b = 0; b < suf_bytes; ++b)
                {
                    uint8_t bval = rec[b];
                    row[base_idx    ] = (bval >> 6u) & 0x3u;
                    row[base_idx + 1] = (bval >> 4u) & 0x3u;
                    row[base_idx + 2] = (bval >> 2u) & 0x3u;
                    row[base_idx + 3] =  bval         & 0x3u;
                    base_idx += 4u;
                }
            }

            // Count — only decode when needed
            if (cnt_ptr)
            {
                uint32_t count = 0;
                for (uint32_t b = 0; b < cnt_bytes; ++b)
                    count |= static_cast<uint32_t>(rec[suf_bytes + b]) << (b * 8u);
                cnt_ptr[i] = count;
            }
        }
    }

    py::dict out;
    out["kmers"]     = std::move(kmers_arr);
    out["key_words"] = std::move(key_words);
    out["n_words"]   = n_words;
    // counts is a zero-element array when drop_count=true — Python layer
    // checks drop_count before accessing raw["counts"].
    out["counts"]    = std::move(counts_arr);
    out["n_unique"]     = r2.nUniqueKmers;
    out["n_total"]      = r2.nTotalKmers;
    out["below_cutoff"] = r2.nBelowCutoffMin;
    out["above_cutoff"] = r2.nAboveCutoffMax;
    out["n_sequences"]  = r1.nSeqences;
    out["stage1_time"]  = r1.time;
    out["stage2_time"]  = r2.time;
    out["kmer_len"]     = kmer_len;
    return out;
}


PYBIND11_MODULE(_core, m)
{
    m.doc() =
        "kmcpy._core — private pybind11 binding for libkmc_core.\n"
        "\n"
        "Stage 1: ramOnlyMode=false — bins in RAM, spill to tmp_path when needed.\n"
        "Stage 2: withoutOutput=true — k-mer counts to RAM, not disk.\n"
        "         strict_mem (default true) — Stage 2 never exceeds max_ram_gb.\n"
        "\n"
        "Output arrays (one pass):\n"
        "  kmers     (n, k)        uint8  — 2-bit codes A=0 C=1 G=2 T=3, all k\n"
        "  key_words (n, n_words)  uint64 — packed key, n_words=ceil(2k/64)\n"
        "  counts    (n,)          uint32 — counts; empty array if drop_count=true\n"
        "  n_words   scalar        uint32 — number of uint64 words per k-mer\n"
        "\n"
        "drop_count=true: counts_arr is allocated with size 0 — no heap allocation,\n"
        "  no fill, no transfer cost.  Saves n*4 bytes at the C++->Python boundary.\n"
        "  Memory saving: 400 MB per 100M k-mers; 3.2 GB for 8 parallel jobs.\n"
        "\n"
        "Space efficiency: n_words is the minimum needed for k bits:\n"
        "  k<=32: n_words=1  k<=64: n_words=2  k<=128: n_words=4  k<=256: n_words=8\n"
        "\n"
        "Key layout: word[0]=MSB (base[0] side), word[n_words-1]=LSB (base[k-1] side)\n"
        "  Right-aligned: base[k-1] at bits 1:0 of word[n_words-1]\n"
        "  Unused high bits of word[0] are always zero.\n"
        "  A=00, C=01, G=10, T=11 (exact KMC 2-bit encoding, zero transformation)";

    m.def(
        "_count_kmers_internal",
        &_count_kmers_internal,
        py::arg("input_files"),
        py::arg("tmp_path"),
        py::arg("kmer_len")      = 25u,
        py::arg("max_ram_gb")    = 12u,
        py::arg("n_threads")     = 0u,
        py::arg("n_readers")     = 0u,
        py::arg("n_splitters")   = 0u,
        py::arg("signature_len") = 9u,
        py::arg("n_bins")        = 0u,
        py::arg("canonical")     = true,
        py::arg("homopolymer")   = false,
        py::arg("strict_mem")    = true,
        py::arg("input_type")    = 0,
        py::arg("cutoff_min")    = 1ULL,
        py::arg("cutoff_max")    = 1000000000ULL,
        py::arg("counter_max")   = 65535u,
        py::arg("drop_count")    = false,
        py::arg("need_bases")    = true,
        "Private — call kmcpy.count_kmers() instead."
    );

    m.attr("__kmc_engine_version__") = KMC::CfgConsts::kmc_ver;
}