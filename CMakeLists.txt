cmake_minimum_required(VERSION 3.4)
project(bucket_mapper CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include(FetchContent)

# seqan3 for the basic skeleton for the mapper
FetchContent_Declare(seqan3
    GIT_REPOSITORY https://github.com/seqan/seqan3.git
    GIT_TAG        cdb6abe5e36b958e86ff1ba94561dac4955d5f14
)
FetchContent_MakeAvailable(seqan3)

# Sharg for argument parsing
FetchContent_Declare(sharg
    GIT_REPOSITORY https://github.com/seqan/sharg-parser.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(sharg)


set(BUILD_FILES
    compressor/bucket_index.hpp)

# build dataset generating tools
add_executable(test main.cpp ${BUILD_FILES})
target_link_libraries(test PUBLIC seqan3::seqan3)

