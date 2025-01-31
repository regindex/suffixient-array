set(SUFFAR_SOURCES suffixient_array_baseline.hpp)

add_library(suffixient_array OBJECT ${SUFFAR_SOURCES})
target_link_libraries(common sdsl divsufsort divsufsort64)