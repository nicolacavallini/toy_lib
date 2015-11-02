add_library(grid hyper_cube.cpp hyper_ball.cpp)

install (TARGETS grid DESTINATION bin)
install (FILES hyper_cube.h DESTINATION include)
install (FILES hyper_ball.h DESTINATION include)
