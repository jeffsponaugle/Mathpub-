// micro-benchmark: ntt_mpn_mul_bonus vs GMP mpn_mul at large sizes, by thread count
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <gmp.h>
#include "ntt.h"
static double now(){return omp_get_wtime();}
int main(int argc,char**argv){
  ntt_init();
  int maxthr = argc>1?atoi(argv[1]):omp_get_max_threads();
  int thr[]={1,2,4,8,10,16,20,32,48,64,96}; int nthr=0; while(nthr<11&&thr[nthr]<=maxthr) nthr++;
  int sizes[]={20,22,24,26}; int nsz = argc>2?atoi(argv[2]):4;
  for(int si=0;si<nsz;si++){
    size_t n=(size_t)1<<sizes[si];
    mp_limb_t*a=malloc(n*8),*b=malloc(n*8),*c=malloc(2*n*8);
    mpn_random(a,n); mpn_random(b,n);
    printf("n=2^%d limbs (%.0f Mbit):",sizes[si],n*64/1e6); fflush(stdout);
    if(sizes[si]<=24){ double t0=now(); mpn_mul(c,a,n,b,n); printf(" gmp=%.2fs",now()-t0); fflush(stdout);}
    mp_limb_t chk=c[n];
    for(int ti=0;ti<nthr;ti++){
      int t=thr[ti]; omp_set_num_threads(t);
      double best=1e9; int reps = sizes[si]<=22?3:1;
      for(int r=0;r<reps;r++){ double t0=now(); ntt_mpn_mul_bonus(c,a,n,1,b,n,1,0,t); double dt=now()-t0; if(dt<best)best=dt; }
      printf(" t%d=%.2fs%s",t,best,(sizes[si]<=24&&c[n]!=chk)?"(MISMATCH)":""); fflush(stdout);
    }
    printf("\n"); free(a);free(b);free(c);
  }
  return 0;
}
