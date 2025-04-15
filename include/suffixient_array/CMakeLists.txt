set(SUFFAR_SOURCES suffixient_array_baseline.hpp suffixient_array_elias_fano.hpp)

add_library(suffixient_array OBJECT ${SUFFAR_SOURCES})
target_link_libraries(common bitvectors sdsl divsufsort divsufsort64)

set(PA_SOURCES prefix_array_binary_search.hpp)

add_library(prefix_array_binary_search OBJECT ${PA_SOURCES})
target_link_libraries(common sdsl divsufsort divsufsort64)