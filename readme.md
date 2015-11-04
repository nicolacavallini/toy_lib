# Toy Lib

Your first libray.

## Installing

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/home/nicola/local/usr/toy_lib ../
make
make install
export TOY_LIB_PATH=/home/nicola/local/usr/toy_lib
```
go to the `tests` directory, to check installation:
```
mkdir build
cd build
cmake ../
ctest
```

## Test suite

### Build and run

Build an run all tests.
```
cd path_to_surce/toy_lib/tests/build
mkdir build
cd build
cmake ../
ctest
```
To run a specific test. Inside the `toy_lib/tests/build`:
```
ctest -V -R feature_dir-test_name
```
The not so smart way to compile and run a test:
```
make feature_dir/test_name.exe
./feature_dir/test_name.exe
```

### Create a new test

- Suppose you want to test the feature `my_feature`. The actual features in `toy_lib` are `grid` and `lac`. 
- Inside `feature_dir` directory create the `test_name_00`. 
- Inside `test_name_00` create an empty `expected.txt` file. 
- Implement your test in the file `test_name_00.cpp`.
- Output the data you want to check in `std::ofstream output_file`, provided in `test.h`, with `output_file << output_I_want_to_check`.
- The desired output will be dumped into `feature_dir/test_name_00/output.txt`, if you ran your test with `ctest -V -R feature_dir-test_name`.
- copy `/output.txt` to `../feature_dir/test_name_00/expected.txt`.
- Control that the test passes:
```
$ ctest -V -R feature_dir-test_name_00

1: Test timeout computed to be: 9.99988e+06
1: ******************************************************
1: ** Compiling and linking test program
1: [100%] Built target feature_dir/test_name_00.exe
1: ** Executing test program.
1: ** Comparing test output with expected one.
1: ******************************************************
1:             Test results
1:             ------------
1: Compiling and linking: 0
1: Runnig program:        0
1: Expected output:       0
1: 
1: Test status: PASSED.
1: ******************************************************
1/1 Test #1: feature_dir-test_name_00 ...............   Passed    0.13 sec

The following tests passed:
	feature_dir-test_name_00

100% tests passed, 0 tests failed out of 1

Total Test time (real) =   0.13 sec
```
