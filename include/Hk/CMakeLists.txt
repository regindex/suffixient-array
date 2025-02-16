set(HK_SOURCES Hk_Ferragina_Venturini.hpp)

add_library(Hk OBJECT ${HK_SOURCES})
target_link_libraries(common bitvectors sdsl divsufsort divsufsort64)