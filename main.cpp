#include "compressor/bucket_index.hpp"


int main() {
    bucket_indexer bi(14, 500);
    bi.read_sequence_file("/home/zhenhao/TDT/data/562.fna");
}