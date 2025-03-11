add_executable(build_test_index build_test_z_fast_index.cpp)
target_include_directories(build_test_index PRIVATE ${CMAKE_SOURCE_DIR}/include/LZ/lz77index/)
target_link_libraries(build_test_index common ctriepp malloc_count ${CMAKE_SOURCE_DIR}/build/lib/liblz77.a divsufsort divsufsort64) 
target_compile_options(build_test_index PRIVATE -Wno-narrowing -w)

add_executable(build_store_index build_store_index.cpp)
target_include_directories(build_store_index PRIVATE ${CMAKE_SOURCE_DIR}/include/LZ/lz77index/)
target_link_libraries(build_store_index common suffixient_array Hk RLZ suffix_array_binary_search malloc_count ${CMAKE_SOURCE_DIR}/build/lib/liblz77.a sdsl divsufsort divsufsort64) 

add_executable(locate locate.cpp)
target_include_directories(locate PRIVATE ${CMAKE_SOURCE_DIR}/include/LZ/lz77index/)
target_link_libraries(locate common suffixient_array Hk RLZ suffix_array_binary_search malloc_count ${CMAKE_SOURCE_DIR}/build/lib/liblz77.a sdsl divsufsort divsufsort64) 

add_executable(mems mems.cpp)
target_include_directories(mems PRIVATE ${CMAKE_SOURCE_DIR}/include/LZ/lz77index/)
target_link_libraries(mems common suffixient_array Hk RLZ suffix_array_binary_search malloc_count ${CMAKE_SOURCE_DIR}/build/lib/liblz77.a sdsl divsufsort divsufsort64) 