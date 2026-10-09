/*
 * symmetry_classes.c -- number of 2n x 2n Takuzu grids up to symmetry, n <= 4,
 * by canonical forms, cross-checked with Burnside's lemma.  Groups: the 8
 * symmetries of the square (D4), and D4 combined with complementation (16).
 *
 * usage: cc -O2 -o symmetry_classes symmetry_classes.c && ./symmetry_classes
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
int N,n,nR; uint32_t R[64]; uint32_t rows[8]; int cnt[8];
uint64_t total, *setD4, *setAll; size_t cap, nD4, nAll; uint64_t fix[16];
static uint64_t enc(int g[8][8]){ uint64_t v=0; for(int i=0;i<N;i++)for(int j=0;j<N;j++) v=(v<<1)|g[i][j]; return v; }
static void tr(int g[8][8], int k, int out[8][8]){ // k: 0-7 dihedral, bit3 complement
  for(int i=0;i<N;i++)for(int j=0;j<N;j++){ int a=i,b=j,t; switch(k&7){
      case 0: break; case 1: t=a; a=b; b=N-1-t; break; case 2: a=N-1-a; b=N-1-b; break; case 3: t=a; a=N-1-b; b=t; break;
      case 4: a=N-1-a; break; case 5: b=N-1-b; break; case 6: t=a; a=b; b=t; break; case 7: t=a; a=N-1-b; b=N-1-t; break; }
    out[i][j]=g[a][b]^((k>>3)&1); } }
static int ins(uint64_t *s, size_t *cnt_, uint64_t v){ size_t h=(v*0x9E3779B97F4A7C15ull)>>20; h%=cap; while(s[h]){ if(s[h]==v+1) return 0; h=(h+1)%cap; } s[h]=v+1; (*cnt_)++; return 1; }
static void leaf(void){ int g[8][8],t[8][8]; for(int i=0;i<N;i++)for(int j=0;j<N;j++) g[i][j]=(rows[i]>>j)&1;
  uint64_t self=enc(g), mD4=~0ull, mAll=~0ull;
  for(int k=0;k<16;k++){ tr(g,k,t); uint64_t v=enc(t); if(v==self) fix[k]++; if(k<8&&v<mD4) mD4=v; if(v<mAll) mAll=v; }
  ins(setD4,&nD4,mD4); ins(setAll,&nAll,mAll); total++; }
static void dfs(int k){ uint32_t mask=(1u<<N)-1;
  if(k==N){ for(int i=0;i<N;i++)for(int j=i+1;j<N;j++) if(rows[i]==rows[j]) return;
    uint32_t c[8]; for(int j=0;j<N;j++){c[j]=0; for(int i=0;i<N;i++) c[j]|=((rows[i]>>j)&1u)<<i;}
    for(int i=0;i<N;i++)for(int j=i+1;j<N;j++) if(c[i]==c[j]) return; leaf(); return; }
  for(int r=0;r<nR;r++){ uint32_t x=R[r]; if(k>=2&&(~(rows[k-2]^rows[k-1])&~(rows[k-1]^x)&mask)) continue;
    int ok=1; for(int j=0;j<N;j++){int cc=cnt[j]+((x>>j)&1); if(cc>n||cc+(N-1-k)<n) ok=0;} if(!ok) continue;
    for(int j=0;j<N;j++) cnt[j]+=(x>>j)&1; rows[k]=x; dfs(k+1); for(int j=0;j<N;j++) cnt[j]-=(x>>j)&1; } }
int main(void){ for(n=1;n<=4;n++){ N=2*n; nR=0; for(uint32_t x=0;x<(1u<<N);x++){ if(__builtin_popcount(x)!=n) continue; int ok=1; for(int i=0;i+2<N;i++){int t=(x>>i)&7; if(t==0||t==7) ok=0;} if(ok) R[nR++]=x; }
    cap=1<<24; setD4=calloc(cap,8); setAll=calloc(cap,8); nD4=nAll=total=0; memset(fix,0,sizeof fix); memset(cnt,0,sizeof cnt); dfs(0);
    uint64_t sD4=0,sAll=0; for(int k=0;k<8;k++) sD4+=fix[k]; for(int k=0;k<16;k++) sAll+=fix[k];
    printf("n=%d grids=%llu  classes under D4: %zu (Burnside %llu)  under D4 x complement: %zu (Burnside %llu)\n", n,(unsigned long long)total,nD4,(unsigned long long)(sD4/8),nAll,(unsigned long long)(sAll/16));
    printf("   fixed points per symmetry (id,r90,r180,r270,flipV,flipH,transpose,antitranspose, then same with complement):"); for(int k=0;k<16;k++) printf(" %llu",(unsigned long long)fix[k]); printf("\n");
    free(setD4); free(setAll); } }
