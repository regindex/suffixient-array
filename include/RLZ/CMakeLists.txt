set(RLZ_SOURCES RLZ_DNA.hpp RLZ_DNA_o2.hpp)

add_library(RLZ OBJECT ${RLZ_SOURCES})
target_link_libraries(common sdsl divsufsort divsufsort64)