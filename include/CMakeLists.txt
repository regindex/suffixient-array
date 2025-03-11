add_subdirectory(common)
target_include_directories(common PUBLIC common)

#add_subdirectory(pfp_iterator)
#target_include_directories(pfp_iterator PUBLIC pfp_iterator)

add_subdirectory(ctriepp)
target_include_directories(ctriepp PUBLIC ctriepp)

add_subdirectory(LZ)
#target_include_directories(LZ PUBLIC LZ)

add_subdirectory(bitvectors)
target_include_directories(bitvectors PUBLIC bitvectors)

add_subdirectory(suffixient_array)
target_include_directories(suffixient_array PUBLIC suffixient_array)

add_subdirectory(Hk)
target_include_directories(Hk PUBLIC Hk)

add_subdirectory(RLZ)
target_include_directories(RLZ PUBLIC RLZ)