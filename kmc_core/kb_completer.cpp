/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.2.4
  Date   : 2024-02-09
*/
#include <algorithm>
#include <numeric>
#include "kb_completer.h"
#include "critical_error_handler.h"
#include <sstream>
// Fork: Required for flat_buf byte accumulator in ProcessBinsFirstStage.
#include <vector>

using namespace std;

extern uint64 total_reads;



//************************************************************************************************************
// CKmerBinCompleter
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Assign queues and monitors
CKmerBinCompleter::CKmerBinCompleter(CKMCParams &Params, CKMCQueues &Queues) 
{
	file_name      = Params.output_file_name;
	kq             = Queues.kq.get();
	bd		       = Queues.bd.get();
	s_mapper	   = Queues.s_mapper.get();
	memory_bins    = Queues.memory_bins.get();

	bbkpq		   = Queues.bbkpq.get();
	use_strict_mem = Params.use_strict_mem;
	kmer_file_name = file_name + ".kmc_suf";
	lut_file_name  = file_name + ".kmc_pre";

	kmer_len       = Params.kmer_len;
	signature_len  = Params.signature_len;

	cutoff_min     = Params.cutoff_min;
	cutoff_max     = (uint32)Params.cutoff_max;
	counter_max    = (uint32)Params.counter_max;
	lut_prefix_len = Params.lut_prefix_len;
	both_strands   = Params.both_strands;
	without_output = Params.without_output;

	kmer_t_size    = Params.KMER_T_size;

	output_type = Params.output_type;		
}


//----------------------------------------------------------------------------------
// Store sorted and compacted bins to the output file (stage first)
void CKmerBinCompleter::ProcessBinsFirstStage()
{
	int32 bin_id = 0;
	uchar *data = nullptr;
	//uint64 data_size = 0;
	list<pair<uint64, uint64>> data_packs;
	uchar *lut = nullptr;
	uint64 lut_size = 0;
	counter_size = 0;
	if (output_type == OutputType::KMC)
	{
		sig_map_size = (1 << (signature_len * 2)) + 1;
		sig_map = new uint32[sig_map_size];
		fill_n(sig_map, sig_map_size, 0);
		lut_pos = 0;
	}
	
	counter_size = calc_counter_size(cutoff_max, counter_max);
	
	if (!without_output)
	{
		if(output_type == OutputType::KMC)
		{
			// Open output file
			out_kmer = fopen(kmer_file_name.c_str(), "wb");
			if (!out_kmer)
			{
				std::ostringstream ostr;
				ostr << "Error: Cannot create " << kmer_file_name;
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}

			out_lut = fopen(lut_file_name.c_str(), "wb");
			if (!out_lut)
			{
				std::ostringstream ostr;
				ostr << "Error: Cannot create " << lut_file_name;
				fclose(out_kmer);
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
		}
		else if (output_type == OutputType::KFF)
		{			
			kff_writer = std::make_unique<CKFFWriter>(file_name + ".kff", both_strands, kmer_len, counter_size, cutoff_min, cutoff_max);
		}
		else
		{
			std::ostringstream ostr;
			ostr << "Error: not implemented, plase contact authors showing this message" << __FILE__ << "\t" << __LINE__;
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
	}
	
	n_recs = 0;

	_n_unique = _n_cutoff_min = _n_cutoff_max = _n_total = 0;
	n_unique  = n_cutoff_min  = n_cutoff_max  = n_total  = 0;

	char s_kmc_pre[] = "KMCP";
	char s_kmc_suf[] = "KMCS";

	if (!without_output)
	{
		if (output_type == OutputType::KMC)
		{
			// Markers at the beginning
			fwrite(s_kmc_pre, 1, 4, out_lut);
			fwrite(s_kmc_suf, 1, 4, out_kmer);
		}
	}

	// Fork: Pre-compute record layout constants for the raw packed accumulator.
	//       kmer_suf_bytes: bytes encoding (kmer_len - lut_prefix_len) suffix bases.
	//       rec_stride:     total bytes per record (suffix bytes + counter bytes).
	//       Used only when without_output == true.
	uint32_t kmer_suf_bytes = 0;
	uint32_t rec_stride     = 0;
	if (without_output && output_type == OutputType::KMC)
	{
		uint32_t suffix_bases = (uint32_t)kmer_len - lut_prefix_len;
		kmer_suf_bytes = (suffix_bases + 3) / 4;
		rec_stride     = kmer_suf_bytes + (uint32_t)counter_size;

		// Store layout metadata into members so GetKmerTable() can pass them
		// to the caller (_core.cpp) for decoding.
		if (!kmer_table_meta_stored)
		{
			kmer_suf_bytes_out     = kmer_suf_bytes;
			counter_size_out       = (uint32_t)counter_size;
			lut_prefix_len_out     = lut_prefix_len;
			kmer_table_meta_stored = true;
		}
	}

	// Process priority queue of ready-to-output bins
	while (!kq->empty())
	{
		// Get the next bin
		if (!kq->pop(bin_id, data, data_packs, lut, lut_size, _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total))
			continue;

		// Decrease memory size allocated by stored bin
		string name;
		uint64 n_rec;
		uint64 n_plus_x_recs;
		uint64 n_super_kmers;
		uint64 raw_size;
		CMemDiskFile *file;

		bd->read(bin_id, file, name, raw_size, n_rec, n_plus_x_recs, n_super_kmers);

		uint64 lut_recs = lut_size / sizeof(uint64);


		if (!without_output)
		{
			if(output_type == OutputType::KMC)
			{ 
				for (auto& e : data_packs)
				{
					// Write bin data to the output file
#ifdef _WIN32 //fwrite bug https://connect.microsoft.com/VisualStudio/feedback/details/755018/fwrite-hangs-with-large-size-count
					uint64 write_offset = e.first;
					uint64 left_to_write = e.second - e.first;
					while (left_to_write)
					{
						uint64 current_to_write = MIN(left_to_write, (4ull << 30) - 1);
						fwrite(data + write_offset, 1, current_to_write, out_kmer);
						write_offset += current_to_write;
						left_to_write -= current_to_write;
					}
#else
					fwrite(data + e.first, 1, e.second - e.first, out_kmer);
#endif
				}
			}
			else if (output_type == OutputType::KFF)
			{
				uint32_t rec_size = (kmer_len + 3) / 4 + counter_size;
				for (auto& e : data_packs)								
					kff_writer->StoreWholeSection(data + e.first, (e.second - e.first) / rec_size);
			}
			else
			{
				std::ostringstream ostr;
				ostr << "Error: not implemented, plase contact authors showing this message" << __FILE__ << "\t" << __LINE__;
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			
		}
		// Fork: When without_output==true, accumulate raw packed bytes into
		//       packed_buf and prefix_buf instead of decoding to ACGT strings.
		//       No ACGT decoding happens here — _core.cpp decodes directly into
		//       numpy buffers, eliminating all intermediate std::string allocations.
		else if (output_type == OutputType::KMC && rec_stride > 0)
		{
			// Flatten data_packs into a contiguous buffer for O(1) record access.
			// data_packs may contain multiple disjoint byte ranges per bin.
			std::vector<uchar> flat_buf;
			{
				uint64_t total_bytes = 0;
				for (auto& e : data_packs)
					total_bytes += (e.second - e.first);
				flat_buf.reserve(total_bytes);
				for (auto& e : data_packs)
					flat_buf.insert(flat_buf.end(), data + e.first, data + e.second);
			}

			const uchar* rec_base = flat_buf.data();
			uint64_t     rec_pos  = 0;

			if (lut_prefix_len == 0)
			{
				// Fork: Special case — lut_prefix_len=0 means no prefix/suffix split.
				// The full k-mer is stored in data_packs as flat records.
				// lut[0] holds the total record count but lut_recs=0, so the normal
				// LUT walk would skip all records. Walk data_packs directly instead.
				// All records get prefix_idx=0 (single slot, full k-mer in rec_base).
				const uint64_t n_recs_bin = flat_buf.size() / rec_stride;
				for (uint64_t r = 0; r < n_recs_bin; ++r)
				{
					if (rec_pos + rec_stride > flat_buf.size())
						break;

					packed_buf.insert(packed_buf.end(),
					                  rec_base + rec_pos,
					                  rec_base + rec_pos + rec_stride);
					prefix_buf.push_back(0);  // prefix_idx=0: full k-mer, no prefix
					rec_pos += rec_stride;
				}
			}
			else
			{
				const uint64_t* ulut = reinterpret_cast<const uint64_t*>(lut);

				// Walk LUT prefix slots and accumulate raw records.
				// Each record = kmer_suf_bytes raw suffix bytes + counter_size count bytes.
				// One prefix_idx stored per record in prefix_buf.
				for (uint64_t prefix_idx = 0; prefix_idx < lut_recs; ++prefix_idx)
				{
					uint64_t delta = ulut[prefix_idx]; // records with this prefix
					for (uint64_t r = 0; r < delta; ++r)
					{
						if (rec_pos + rec_stride > flat_buf.size())
							break; // guard: malformed bin, should not happen

						// Fork: Append raw bytes as-is — no decoding.
						packed_buf.insert(packed_buf.end(),
						                  rec_base + rec_pos,
						                  rec_base + rec_pos + rec_stride);

						// Fork: Store prefix index for this record.
						prefix_buf.push_back(prefix_idx);

						rec_pos += rec_stride;
					}
				}
			}
		}

		memory_bins->free(bin_id, CMemoryBins::mba_suffix);

		if (!without_output)
		{
			if (output_type == OutputType::KMC)
			{
				uint64* ulut = (uint64*)lut;
				for (uint64 i = 0; i < lut_recs; ++i)
				{
					uint64 x = ulut[i];
					ulut[i] = n_recs;
					n_recs += x;
				}
				fwrite(lut, lut_recs, sizeof(uint64), out_lut);
			}
		}
		//fwrite(&n_rec, 1, sizeof(uint64), out_lut);
		memory_bins->free(bin_id, CMemoryBins::mba_lut);

		n_unique	 += _n_unique;
		n_cutoff_min += _n_cutoff_min;
		n_cutoff_max += _n_cutoff_max;
		n_total      += _n_total;
		
		if (output_type == OutputType::KMC)
		{
			for (uint32 i = 0; i < sig_map_size; ++i)
			{
				if (s_mapper->get_bin_id(i) == bin_id)
				{
					sig_map[i] = lut_pos;
				}
			}
			++lut_pos;
		}
	}		
}

//----------------------------------------------------------------------------------
// Store sorted and compacted bins to the output file (stage second)
void CKmerBinCompleter::ProcessBinsSecondStage()
{
	char s_kmc_pre[] = "KMCP";
	char s_kmc_suf[] = "KMCS";
	if (use_strict_mem)
	{
		int32 bin_id;
		uchar *data = nullptr;
		uint64 data_size = 0;
		uchar *lut = nullptr;
		uint64 lut_size = 0;		
		bool last_in_bin = false;
		while (bbkpq->pop(bin_id, data, data_size, lut, lut_size, _n_unique, _n_cutoff_min, _n_cutoff_max, _n_total, last_in_bin))
		{
			if (data_size)
			{
				if(!without_output)
					fwrite(data, 1, data_size, out_kmer);				
				sm_pmm_merger_suff->free(data);
			}
			if (lut_size)
			{
				uint64 lut_recs = lut_size / sizeof(uint64);
				uint64* ulut = (uint64*)lut;
				for (uint64 i = 0; i < lut_recs; ++i)
				{
					uint64 x = ulut[i];
					ulut[i] = n_recs;
					n_recs += x;
				}
				if(!without_output)
					fwrite(lut, lut_recs, sizeof(uint64), out_lut);				
				sm_pmm_merger_lut->free(lut);
			}
			if(last_in_bin)
			{
				n_unique += _n_unique;
				n_cutoff_min += _n_cutoff_min;
				n_cutoff_max += _n_cutoff_max;
				n_total += _n_total;
				for (uint32 i = 0; i < sig_map_size; ++i)
				{
					if (s_mapper->get_bin_id(i) == bin_id)
					{
						sig_map[i] = lut_pos;
					}
				}
				++lut_pos;
			}
		}
	}

	if (!without_output)
	{
		if(output_type == OutputType::KMC)
		{
			// Marker at the end
			fwrite(s_kmc_suf, 1, 4, out_kmer);
			fclose(out_kmer);

			fwrite(&n_recs, 1, sizeof(uint64), out_lut);

			//store signature mapping 
			fwrite(sig_map, sizeof(uint32), sig_map_size, out_lut);

			// Store header
			uint32 offset = 0;

			store_uint(out_lut, kmer_len, 4);				offset += 4;
			store_uint(out_lut, (uint32)0, 4);				offset += 4;	// mode: 0 (counting), 1 (Quake-compatibile counting) which is now not supported
			store_uint(out_lut, counter_size, 4);			offset += 4;
			store_uint(out_lut, lut_prefix_len, 4);			offset += 4;
			store_uint(out_lut, signature_len, 4);			offset += 4;
			store_uint(out_lut, cutoff_min, 4);				offset += 4;
			store_uint(out_lut, cutoff_max, 4);				offset += 4;
			store_uint(out_lut, n_unique - n_cutoff_min - n_cutoff_max, 8);		offset += 8;

			store_uint(out_lut, both_strands ? 0 : 1, 1);			offset++;

			// Space for future use
			for (int32 i = 0; i < 27; ++i)
			{
				store_uint(out_lut, 0, 1);
				offset++;
			}

			store_uint(out_lut, 0x200, 4);
			offset += 4;

			store_uint(out_lut, offset, 4);

			// Marker at the end
			fwrite(s_kmc_pre, 1, 4, out_lut);
			fclose(out_lut);
		}
	}

	if (output_type == OutputType::KMC)	
		delete[] sig_map;
}

//----------------------------------------------------------------------------------
// Return statistics
void CKmerBinCompleter::GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total)
{
	_n_unique	  = n_unique;
	_n_cutoff_min = n_cutoff_min;
	_n_cutoff_max = n_cutoff_max;
	_n_total      = n_total;
}

//----------------------------------------------------------------------------------
// Store single unsigned integer in LSB fashion
bool CKmerBinCompleter::store_uint(FILE *out, uint64 x, uint32 size)
{
	for(uint32 i = 0; i < size; ++i)
		putc((x >> (i * 8)) & 0xFF, out);

	return true;
}

//----------------------------------------------------------------------------------
//Init memory pools for 2nd stage
void CKmerBinCompleter::InitStage2(CKMCParams& /*Params*/, CKMCQueues& Queues)
{
	sm_pmm_merger_lut = Queues.sm_pmm_merger_lut.get();
	sm_pmm_merger_suff = Queues.sm_pmm_merger_suff.get();
}


//************************************************************************************************************
// CWKmerBinCompleter
//************************************************************************************************************

//----------------------------------------------------------------------------------
// Constructor
CWKmerBinCompleter::CWKmerBinCompleter(CKMCParams &Params, CKMCQueues &Queues)
{
	kbc = std::make_unique<CKmerBinCompleter>(Params, Queues);
}

void CWKmerBinCompleter::InitStage2(CKMCParams& Params, CKMCQueues& Queues)
{
	kbc->InitStage2(Params, Queues);
}

//----------------------------------------------------------------------------------
// Execution
void CWKmerBinCompleter::operator()(bool first_stage)
{
	if(first_stage)
		kbc->ProcessBinsFirstStage();
	else
		kbc->ProcessBinsSecondStage();
}

//----------------------------------------------------------------------------------
// Return statistics
void CWKmerBinCompleter::GetTotal(uint64 &_n_unique, uint64 &_n_cutoff_min, uint64 &_n_cutoff_max, uint64 &_n_total)
{
	if(kbc)
		kbc->GetTotal(_n_unique, _n_cutoff_min, _n_cutoff_max, _n_total);
}

// ***** EOF