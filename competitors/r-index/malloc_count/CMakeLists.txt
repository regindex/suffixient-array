set(MALLOCOUNT_SOURCES malloc_count.c)

add_library(malloc_count OBJECT ${MALLOCOUNT_SOURCES})
target_link_libraries(dl)
