rm *.o
make -f Makefile_MPFR.gcc
mv SweeD SweeD-MPFR 
rm *.o;

make -f Makefile.gcc && rm *.o;

make -f ./Makefile_MPFR.PTHREADS.gcc
mv SweeD-P SweeD-MPFR-P
rm *.o;


make -f ./Makefile.PTHREADS.gcc
rm *.o;
