export PSTL_USAGE_WARNINGS=1
export ONEDPL_USE_DPCPP_BACKEND=1

dpcpp -std=c++17 -O3 src/propagate-tor-test_sycl.cpp -o test-dpl.exe -Dntrks=8192 -Dnevts=100 -DNITER=5 -Dbsize=1 -Dnlayer=20


