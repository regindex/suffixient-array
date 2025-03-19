add_subdirectory(common_)
target_include_directories(common_ PUBLIC common_)

add_subdirectory(pfp_iterator)
target_include_directories(pfp_iterator PUBLIC pfp_iterator)

add_subdirectory(ctriepp)
target_include_directories(ctriepp PUBLIC ctriepp)

add_subdirectory(LZ)

add_subdirectory(bitvectors)
target_include_directories(bitvectors PUBLIC bitvectors)

add_subdirectory(suffixient_array)
target_include_directories(suffixient_array PUBLIC suffixient_array)

add_subdirectory(RLZ)
target_include_directories(RLZ PUBLIC RLZ)