add_executable(pscan pscan.cpp utils.c xerrors.c)
target_link_libraries(pscan pthread) 