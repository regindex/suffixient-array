set(SUFFAR_SOURCES suffixient_array_baseline.hpp suffixient_array_elias_fano.hpp)

add_library(suffixient_array OBJECT ${SUFFAR_SOURCES})
target_link_libraries(common bitvectors sdsl divsufsort divsufsort64)

set(SA_SOURCES suffix_array_binary_search.hpp)

add_library(suffix_array_binary_search OBJECT ${SA_SOURCES})
target_link_libraries(common sdsl divsufsort divsufsort64)