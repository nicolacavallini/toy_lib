find_package(BLAS)

if(BLAS_FOUND)
    add_definitions(-DUSEBLAS)
endif(BLAS_FOUND)


add_library(lac matrix.cpp)

install (TARGETS lac DESTINATION bin)
install (FILES matrix.h DESTINATION include)
