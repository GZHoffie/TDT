#ifndef BUCKET_COMPRESSOR_BUCKET_INDEX_HPP
#define BUCKET_COMPRESSOR_BUCKET_INDEX_HPP

#include <string>
#include <vector>
#include <map>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
 

class bucket_indexer {
private:
    unsigned int k; // length of the minimizer
    unsigned int l; // length of window to take minimizer
    uint64_t seed = 0x8F3F73B5CF1C9ADE; // seed used by seqan3 to find minimizer

    typedef std::pair<uint64_t, uint64_t> node_t;

    std::map<node_t, unsigned int> graph;


    /**
     * @brief turns a sequence into a vector of minimizers.
     */

    std::vector<uint64_t> _to_minimizers(seqan3::dna4_vector& sequence) {
        auto minimisers = sequence | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{k}}, seqan3::window_size{l});
        auto hash_values = minimisers | std::views::transform([&](uint64_t i){ return i ^ seed; });
        std::vector<uint64_t> res(hash_values.begin(), hash_values.end());
        return res;
    }

public:



};



#endif