set(PFPIT_SOURCES pfp_iterator.hpp pfp.hpp common.hpp)

add_library(pfp_iterator OBJECT ${PFPIT_SOURCES})
target_link_libraries(sdsl divsufsort divsufsort64)
#target_link_libraries(pfp common sdsl divsufsort divsufsort64)