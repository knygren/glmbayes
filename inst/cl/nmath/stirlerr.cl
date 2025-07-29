// stirlerr.cl - OpenCL Adaptation of stirlerr.c
//@provides: stirlerr
//@depends: nmath, stirlerr_small, stirlerr_large
//@includes: nmath


inline double stirlerr(double n) {

     
        if (n <= 5.25)         return stirlerr_small(n);
        else return stirlerr_large(n);

} 














