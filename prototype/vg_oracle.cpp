// Reference output from gbwtgraph itself, for testing the Python prototype.
// Reads lines from stdin:
//   key K KMER        -> forward key, forward hash, reverse-complement key, reverse-complement hash
//   min K W SEQUENCE  -> one line per minimizer: offset, strand (+/-), key, hash
// (gbwtgraph reports a reverse-strand minimizer at the offset of its last base.)
// Build (link against the vg build tree):
//   g++ -O2 -std=c++17 -I$VG/include vg_oracle.cpp -o vg_oracle -L$VG/lib \
//       -lgbwtgraph -lgbwt -lhandlegraph -lsdsl -ldivsufsort -ldivsufsort64 -fopenmp
#include <gbwtgraph/minimizer.h>

#include <iostream>
#include <sstream>
#include <string>

using gbwtgraph::DefaultMinimizerIndex;
using gbwtgraph::Key64;

int main()
{
  std::string line;
  while(std::getline(std::cin, line))
  {
    std::istringstream fields(line);
    std::string mode;
    size_t k = 0;
    fields >> mode >> k;
    if(mode == "key")
    {
      std::string kmer;
      fields >> kmer;
      Key64 forward = Key64::encode(kmer);
      Key64 reverse = forward.reverse_complement(k);
      std::cout << forward.get_key() << "\t" << forward.hash() << "\t"
                << reverse.get_key() << "\t" << reverse.hash() << "\n";
    }
    else if(mode == "min")
    {
      size_t w = 0;
      std::string sequence;
      fields >> w >> sequence;
      DefaultMinimizerIndex index(k, w, 1);
      for(auto& minimizer : index.minimizers(sequence))
      {
        std::cout << minimizer.offset << "\t" << (minimizer.is_reverse ? "-" : "+") << "\t"
                  << minimizer.key.get_key() << "\t" << minimizer.hash << "\n";
      }
      std::cout << "end\n";
    }
  }
  return 0;
}
