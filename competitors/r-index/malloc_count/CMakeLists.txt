set(MALLOCOUNT_SOURCES malloc_count.c)

add_library(malloc_count OBJECT ${MALLOCOUNT_SOURCES})
find_library(DL_LIBRARY NAMES dl)
target_link_libraries(malloc_count PRIVATE ${DL_LIBRARY})

#target_link_libraries(dl)
