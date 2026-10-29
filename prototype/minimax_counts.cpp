// The count table minimaxmer selection needs, for every canonical k-mer in the graph.
// Minimaxmers come from Yoshi Suzuki's himapper_rp (https://github.com/yoshihikosuzuki/himapper_rp).
//
//   genome_count  number of distinct graph positions of the k-mer over all haplotypes: the count that
//                 gbwtgraph::frequent_kmers() thresholds on (values.second of build_kmer_indexes()), and the
//                 graph analogue of Yoshi's occurrences in the genome.
//   reads_count   count of the k-mer in a KFF file (KMC on the sample's reads), 0 if absent.
//
// Both saturate at 65535. A k-mer is canonical as in gbwtgraph: the smaller 2-bit key of its two strands.
//
// Usage: minimax_counts -g GRAPH.gbz -o PREFIX [-f READS.kff] [-k K] [-c COVERAGE] [-t THREADS] [--dump-text]
//   -k is required without -f; with -f it must match the KFF k if given.
// Writes (numbers only, never sequences):
//   PREFIX.mmc            counts file: "MMXCNT01", uint64 k, uint64 n, n uint64 keys (sorted Key64::get_key()),
//                         n uint32 (genome_count << 16 | reads_count); little-endian
//   PREFIX.summary.tsv    headline numbers first, then the rest
//   PREFIX.genome_hist.tsv, PREFIX.reads_hist.tsv, PREFIX.joint.tsv
//   PREFIX.dump.tsv       with --dump-text (tests only): key, genome_count, reads_count per k-mer
//
// Build (link against the vg build tree; kff.cpp is vg's ParallelKFFReader):
//   g++ -O3 -std=c++17 -fopenmp -I$VG/include -I$VG/src minimax_counts.cpp $VG/src/kff.cpp -o minimax_counts \
//       -L$VG/lib -lgbwtgraph -lgbwt -l:libhandlegraph.a -lsdsl -ldivsufsort -ldivsufsort64 -lkff -lcrypto -lzstd
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/index.h>
#include <gbwt/utils.h>
#include "kff.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using gbwtgraph::Key64;
using index_type = gbwtgraph::KmerIndex<Key64>;

constexpr std::uint64_t COUNT_MAX = 65535;
constexpr char MAGIC[8] = { 'M', 'M', 'X', 'C', 'N', 'T', '0', '1' };
constexpr size_t KFF_BLOCK_SIZE = 1 << 20;
constexpr std::uint64_t PEAK_MIN_READS = 5;  // the reads_count peak ignores the low (error / absent) end

struct Config
{
  std::string graph, kff, prefix;
  size_t k = 0, threads = 0;
  double coverage = 0.0;
  bool dump_text = false;
};

[[noreturn]] void fail(const std::string& message)
{
  std::cerr << "minimax_counts: " << message << std::endl;
  std::exit(1);
}

Config parse(int argc, char** argv)
{
  Config config;
  for(int i = 1; i < argc; i++)
  {
    std::string arg = argv[i];
    auto value = [&]() -> std::string { if(i + 1 >= argc) { fail("missing value for " + arg); } return argv[++i]; };
    if(arg == "-g") { config.graph = value(); }
    else if(arg == "-f") { config.kff = value(); }
    else if(arg == "-o") { config.prefix = value(); }
    else if(arg == "-k") { config.k = std::stoul(value()); }
    else if(arg == "-c") { config.coverage = std::stod(value()); }
    else if(arg == "-t") { config.threads = std::stoul(value()); }
    else if(arg == "--dump-text") { config.dump_text = true; }
    else { fail("unknown argument " + arg); }
  }
  if(config.graph.empty() || config.prefix.empty()) { fail("usage: -g GRAPH.gbz -o PREFIX [-f READS.kff] [-k K] [-c COVERAGE] [-t THREADS] [--dump-text]"); }
  if(config.kff.empty() && config.k == 0) { fail("-k is required without -f"); }
  return config;
}

// Initial hash table size as in vg minimizer: half the length of the reference and non-alt generic paths
// (or of all nodes if there are no such paths) per index, so memory matches a weighted index build.
size_t hash_table_size(const gbwtgraph::GBWTGraph& graph)
{
  size_t genome_size = 0;
  graph.for_each_path_of_sense({ handlegraph::PathSense::REFERENCE, handlegraph::PathSense::GENERIC },
    [&](const handlegraph::path_handle_t& path)
  {
    if(graph.get_path_name(path).compare(0, 5, "_alt_") == 0) { return; }
    graph.for_each_step_in_path(path, [&](const handlegraph::step_handle_t& step)
    {
      genome_size += graph.get_length(graph.get_handle_of_step(step));
    });
  });
  if(genome_size == 0)
  {
    graph.for_each_handle([&](const handlegraph::handle_t& handle) { genome_size += graph.get_length(handle); });
  }
  std::cerr << "Genome size for the hash table estimate: " << genome_size << std::endl;
  return index_type::minimum_size(genome_size / 2);
}

std::uint64_t canonical(std::uint64_t key, size_t k)
{
  return std::min(key, Key64(key).reverse_complement(k).get_key());
}

void write_hist(const std::string& path, const std::string& column, const std::map<std::uint64_t, std::uint64_t>& hist)
{
  std::ofstream out(path);
  out << column << "\tkmers\n";
  for(auto& [value, count] : hist) { out << value << "\t" << count << "\n"; }
}

// Bin label for the joint table.
std::string bin(std::uint64_t value, const std::vector<std::uint64_t>& upper, const std::vector<std::string>& labels)
{
  for(size_t i = 0; i < upper.size(); i++) { if(value <= upper[i]) { return labels[i]; } }
  return labels.back();
}

int main(int argc, char** argv)
{
  Config config = parse(argc, argv);
  if(config.threads > 0) { omp_set_num_threads(config.threads); }
  double start = gbwt::readTimer();

  // Open the KFF first so that a wrong k fails before the expensive graph pass.
  std::unique_ptr<vg::ParallelKFFReader> reader;
  if(!config.kff.empty())
  {
    reader.reset(new vg::ParallelKFFReader(config.kff));
    if(config.k != 0 && config.k != reader->k) { fail("-k " + std::to_string(config.k) + " but the KFF has k=" + std::to_string(reader->k)); }
    config.k = reader->k;
    std::cerr << "KFF: k=" << reader->k << ", count bytes " << reader->data_bytes << std::endl;
  }
  if(config.k > Key64::KMER_MAX_LENGTH) { fail("k must be at most " + std::to_string(Key64::KMER_MAX_LENGTH)); }

  gbwtgraph::GBZ gbz;
  sdsl::simple_sds::load_from(gbz, config.graph);
  std::cerr << "Loaded the graph in " << (gbwt::readTimer() - start) << " seconds" << std::endl;

  // Genome counts: the distinct graph positions of every canonical k-mer.
  double checkpoint = gbwt::readTimer();
  size_t table_size = hash_table_size(gbz.graph);
  std::array<index_type, 4> indexes
  {
    index_type(table_size, 0), index_type(table_size, 0), index_type(table_size, 0), index_type(table_size, 0)
  };
  gbwtgraph::build_kmer_indexes(gbz.graph, indexes, config.k);
  size_t total = 0;
  for(auto& index : indexes) { total += index.size(); }
  std::cerr << "Counted " << total << " distinct k-mers in the graph in " << (gbwt::readTimer() - checkpoint) << " seconds" << std::endl;

  std::vector<std::pair<std::uint64_t, std::uint32_t>> table;
  table.reserve(total);
  for(auto& index : indexes)
  {
    index.for_each_kmer([&](Key64 key, gbwtgraph::KmerEncoding::multi_value_type values)
    {
      table.emplace_back(key.get_key(), std::min<std::uint64_t>(values.second, COUNT_MAX));
    });
    index = index_type(index_type::INITIAL_CAPACITY, 0);  // free the positions as soon as possible
  }
  gbwt::parallelQuickSort(table.begin(), table.end());
  for(size_t i = 1; i < table.size(); i++)
  {
    if(table[i].first == table[i - 1].first) { fail("a k-mer occurs in two middle-base indexes"); }
  }

  // Reads counts: stream the KFF and add counts for k-mers in the graph; both strands map to the canonical key.
  std::vector<std::uint32_t> reads(table.size(), 0);
  std::uint64_t kff_kmers = 0, kff_in_graph = 0, kff_occurrences = 0, kff_occurrences_in_graph = 0;
  if(reader)
  {
    checkpoint = gbwt::readTimer();
    #pragma omp parallel reduction(+:kff_kmers, kff_in_graph, kff_occurrences, kff_occurrences_in_graph)
    {
      while(true)
      {
        std::vector<std::pair<vg::ParallelKFFReader::kmer_type, size_t>> block = reader->read(KFF_BLOCK_SIZE);
        for(auto& [kmer, count] : block)
        {
          std::uint64_t key = canonical(kmer, config.k);
          auto iter = std::lower_bound(table.begin(), table.end(), std::make_pair(key, std::uint32_t(0)));
          kff_kmers++; kff_occurrences += count;
          if(iter != table.end() && iter->first == key)
          {
            __atomic_fetch_add(&reads[iter - table.begin()], static_cast<std::uint32_t>(std::min<size_t>(count, COUNT_MAX)), __ATOMIC_RELAXED);
            kff_in_graph++; kff_occurrences_in_graph += count;
          }
        }
        if(block.size() < KFF_BLOCK_SIZE) { break; }
      }
    }
    std::cerr << "Read " << kff_kmers << " KFF k-mers in " << (gbwt::readTimer() - checkpoint) << " seconds" << std::endl;
  }

  // Counts file.
  std::string counts_path = config.prefix + ".mmc";
  {
    std::ofstream out(counts_path, std::ios::binary);
    std::uint64_t k = config.k, n = table.size();
    out.write(MAGIC, sizeof(MAGIC));
    out.write(reinterpret_cast<const char*>(&k), sizeof(k));
    out.write(reinterpret_cast<const char*>(&n), sizeof(n));
    for(auto& entry : table) { out.write(reinterpret_cast<const char*>(&entry.first), sizeof(entry.first)); }
    for(size_t i = 0; i < table.size(); i++)
    {
      std::uint32_t packed = (table[i].second << 16) | std::min<std::uint32_t>(reads[i], COUNT_MAX);
      out.write(reinterpret_cast<const char*>(&packed), sizeof(packed));
    }
    if(!out) { fail("cannot write " + counts_path); }
  }
  std::uint64_t counts_bytes = std::filesystem::file_size(counts_path);

  // Distributions.
  std::map<std::uint64_t, std::uint64_t> genome_hist, reads_hist, reads_hist_unique;
  std::map<std::pair<std::string, std::string>, std::uint64_t> joint;
  const std::vector<std::uint64_t> genome_upper { 1, 2, 10, 100 };
  const std::vector<std::string> genome_labels { "1", "2", "3-10", "11-100", ">100" };
  const std::vector<std::uint64_t> reads_upper { 0, 1, 9, 19, 49, 99, 199, 999 };
  const std::vector<std::string> reads_labels { "0", "1", "2-9", "10-19", "20-49", "50-99", "100-199", "200-999", ">=1000" };
  std::uint64_t genome_unique = 0, saturated_genome = 0, saturated_reads = 0, reads_zero = 0;
  for(size_t i = 0; i < table.size(); i++)
  {
    std::uint64_t g = table[i].second, r = std::min<std::uint32_t>(reads[i], COUNT_MAX);
    genome_hist[g]++; reads_hist[r]++;
    if(g == 1) { genome_unique++; reads_hist_unique[r]++; }
    if(g == COUNT_MAX) { saturated_genome++; }
    if(r == COUNT_MAX) { saturated_reads++; }
    if(r == 0) { reads_zero++; }
    joint[{ bin(g, genome_upper, genome_labels), bin(r, reads_upper, reads_labels) }]++;
  }
  std::uint64_t peak = 0, peak_kmers = 0;
  for(auto& [r, count] : reads_hist_unique)
  {
    if(r >= PEAK_MIN_READS && r < COUNT_MAX && count > peak_kmers) { peak = r; peak_kmers = count; }
  }
  double n = static_cast<double>(std::max<size_t>(table.size(), 1));

  write_hist(config.prefix + ".genome_hist.tsv", "genome_count", genome_hist);
  write_hist(config.prefix + ".reads_hist.tsv", "reads_count", reads_hist);
  {
    std::ofstream out(config.prefix + ".joint.tsv");
    out << "genome_count\treads_count\tkmers\n";
    for(auto& g : genome_labels)
    {
      for(auto& r : reads_labels)
      {
        auto iter = joint.find({ g, r });
        out << g << "\t" << r << "\t" << (iter == joint.end() ? 0 : iter->second) << "\n";
      }
    }
  }
  {
    std::ofstream out(config.prefix + ".summary.tsv");
    out << "metric\tvalue\n";
    // Headline numbers first.
    out << "distinct_graph_kmers\t" << table.size() << "\n";
    out << "genome_count_1_fraction\t" << genome_unique / n << "\n";
    if(reader)
    {
      out << "reads_count_peak (genome_count 1, reads_count >= " << PEAK_MIN_READS << ")\t" << peak << "\n";
      if(config.coverage > 0) { out << "reads_count_peak / coverage " << config.coverage << "\t" << peak / config.coverage << "\n"; }
    }
    out << "counts_file_bytes\t" << counts_bytes << "\n";
    // The rest.
    out << "k\t" << config.k << "\n";
    out << "genome_count_saturated\t" << saturated_genome << "\n";
    if(reader)
    {
      out << "reads_count_0_fraction\t" << reads_zero / n << "\n";
      out << "reads_count_saturated\t" << saturated_reads << "\n";
      out << "kff_kmers\t" << kff_kmers << "\n";
      out << "kff_kmers_in_graph\t" << kff_in_graph << "\n";
      out << "kff_kmers_not_in_graph\t" << (kff_kmers - kff_in_graph) << "\n";
      out << "kff_occurrences\t" << kff_occurrences << "\n";
      out << "kff_occurrences_in_graph\t" << kff_occurrences_in_graph << "\n";
    }
    out << "seconds\t" << (gbwt::readTimer() - start) << "\n";
    out << "peak_rss_gib\t" << gbwt::inGigabytes(gbwt::memoryUsage()) << "\n";
  }
  if(config.dump_text)
  {
    std::ofstream out(config.prefix + ".dump.tsv");
    out << "key\tgenome_count\treads_count\n";
    for(size_t i = 0; i < table.size(); i++)
    {
      out << table[i].first << "\t" << table[i].second << "\t" << std::min<std::uint32_t>(reads[i], COUNT_MAX) << "\n";
    }
  }
  std::cerr << "Done in " << (gbwt::readTimer() - start) << " seconds" << std::endl;
  return 0;
}
