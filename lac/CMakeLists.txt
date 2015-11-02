add_library(lac matrix.cpp)

install (TARGETS lac DESTINATION bin)
install (FILES matrix.h DESTINATION include)
