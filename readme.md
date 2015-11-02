# Toy Lib

Your first libray.

```
cmake -DCMAKE_INSTALL_PREFIX=/home/nicola/local/usr/toy_lib ../
make
make install
export TOY_LIB_PATH=/home/nicola/local/usr/toy_lib
```
go to the `tests` directory:
```
cmake ../
ctest
```