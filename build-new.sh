#export GCC_HOME=/home/astrel/install/spack-0.16.3/opt/spack/linux-ubuntu20.04-skylake/gcc-9.3.0/gcc-12.1.0-ocudpow4q7xwavge4c4lzyl5z7kiz5em

export GCC_HOME=/home/astrel/install/spack-0.16.3/opt/spack/linux-ubuntu20.04-skylake/gcc-9.3.0/gcc-11.2.0-wby54stmjhxbuxrgsxwjz3zzr4edo6o4

SRC=./src/propagate-tor-test_cuda_hybrid_cpp2x_v2.cpp

nvc++ -O2 -std=c++23 --gcc-toolchain=$GCC_HOME -stdpar=gpu -gpu=cc86 -gpu=managed -gpu=fma -gpu=fastmath -gpu=autocollapse -gpu=loadcache:L1 -gpu=unroll $SRC  -o ./propagate_tor_stdexec -Dntrks=8192 -Dnevts=100 -DNITER=100 -Dbsize=32 -Dnlayer=20 # -Dinclude_data

nvc++ -O2 -std=c++23 --gcc-toolchain=$GCC_HOME -stdpar=gpu -gpu=cc86 -gpu=managed -gpu=fma -gpu=fastmath -gpu=autocollapse -gpu=loadcache:L1 -gpu=unroll $SRC  -o ./propagate_tor_cuda -Dntrks=8192 -Dnevts=100 -DNITER=100 -Dbsize=32 -Dnlayer=20  -Dcuda_launcher #-Dinclude_data 


#nvc++ -O2 -std=c++23 --gcc-toolchain=$GCC_HOME -stdpar=gpu -gpu=cc86 -gpu=managed -gpu=fma -gpu=autocollapse -gpu=loadcache:L1 -gpu=unroll ./src/propagate-tor-test_cuda_hybrid_mdspan.cpp  -o ./propagate_nvcpp_cuda_nofma -Dntrks=8192 -Dnevts=100 -DNITER=5 -Dbsize=32 -Dnlayer=20 -Dinclude_data

#nvc++ -O2 -std=c++20 --gcc-toolchain=$GCC_HOME -stdpar=gpu -gpu=cc86 -gpu=managed -gpu=fma -gpu=fastmath -gpu=autocollapse -gpu=loadcache:L1 -gpu=unroll ./src/propagate-tor-test_cuda_hybrid.cpp  -o ./propagate_nvcpp_cuda_nvcpp -Dntrks=8192 -Dnevts=100 -DNITER=5 -Dbsize=32 -Dnlayer=20  #-Dstdpar_launcher #-Dinclude_data

#nvc++ -O2 -std=c++20 --gcc-toolchain=$GCC_HOME -stdpar=multicore  ./src/propagate-tor-test_cuda_hybrid.cpp  -o ./propagate_nvcpp_x86_nvcpp -Dntrks=8192 -Dnevts=100 -DNITER=5 -Dbsize=32 -Dnlayer=20 #-Dinclude_data

