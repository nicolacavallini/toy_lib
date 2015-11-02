add_library(grid hyper_cube.cpp)

install (TARGETS grid DESTINATION bin)
install (FILES hyper_cube.h DESTINATION include)
