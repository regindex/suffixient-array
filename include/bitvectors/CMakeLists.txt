set(BITVECTORS_SOURCES succinct_bitvector.hpp elias_fano_bitvector.hpp)

add_library(bitvectors OBJECT ${BITVECTORS_SOURCES})
target_link_libraries(common_ sdsl divsufsort divsufsort64)