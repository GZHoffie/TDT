#ifndef BUCKET_COMPRESSOR_BUCKET_INDEX_HPP
#define BUCKET_COMPRESSOR_BUCKET_INDEX_HPP

#include <string>
#include <vector>
#include <map>
#include <filesystem>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/io/sequence_file/input.hpp>
 

class bucket_indexer {
private:
    uint8_t k; // length of the minimizer
    unsigned int l; // length of window to take minimizer
    uint64_t seed = 0x8F3F73B5CF1C9ADE; // seed used by seqan3 to find minimizer

    typedef std::pair<uint64_t, uint64_t> node_t;

    std::map<node_t, unsigned int> graph;


    /**
     * @brief turns a sequence into a vector of minimizers.
     */

    std::vector<uint64_t> _to_minimizers(const seqan3::dna4_vector& sequence) {
        auto minimisers = sequence | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{k}}, seqan3::window_size{l});
        auto hash_values = minimisers | std::views::transform([&](uint64_t i){ return i ^ seed; });
        //seqan3::debug_stream << hash_values;
        std::vector<uint64_t> res;
        for (auto hash: hash_values) {
            res.push_back(hash);
        }
        return res;
    }

public:
    bucket_indexer(uint8_t minimizer_length, unsigned int window_size) {
        k = minimizer_length;
        l = window_size;
    }

    void read(const seqan3::dna4_vector& sequence) {
        auto minimizers = _to_minimizers(sequence);
        for (unsigned int i = 0; i < minimizers.size() - 1; i++) {
            graph[{minimizers[i], minimizers[i+1]}] += 1;
        }
        seqan3::debug_stream << minimizers.size() << "\n";
        seqan3::debug_stream << graph.size() << "\n";
        //seqan3::debug_stream << minimizers << "\n";
    }

    void read_sequence_file(const std::filesystem::path& fasta_file) {
        seqan3::sequence_file_input fin{fasta_file};

        unsigned int index = 0;

        for (auto & record : fin) {
            seqan3::debug_stream << record.sequence().size() << "\n";

            std::vector<seqan3::dna4> seq(record.sequence().begin(), record.sequence().end());
            read(seq);





            
            index++;
            if (index >= 2) {
                break;
            }
        }
    }
};



#endif