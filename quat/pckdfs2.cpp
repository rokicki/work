#include <string.h>
struct puz {
   char corner[20] ;
   char center[12] ;
   bool operator<(const puz &p) const {
      return memcmp(this, &p, sizeof(puz))<0 ;
   }
   bool operator==(const puz &p) const {
      return memcmp(this, &p, sizeof(puz))==0 ;
   }
   bool operator!=(const puz &p) const {
      return memcmp(this, &p, sizeof(puz))!=0 ;
   }
} ;
void init(puz &p) {
   p = puz({
         {0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57},
         {0,1,2,3,4,5,6,7,8,9,10,11},
       }) ;
}
int inc(int a, int b, int m) {
   return a - a % m + (a + b) % m ;
}
void slowmove(const puz &src, puz &dst, int m, int n) {
   puz t = src ;
   dst = src ;
   while (n-- > 0) {
      switch (m) {
case 0:
         dst.corner[5] = t.corner[15] ;
         dst.corner[6] = inc(t.corner[16],2,3) ;
         dst.corner[8] = t.corner[13] ;
         dst.corner[9] = t.corner[5] ;
         dst.corner[10] = t.corner[8] ;
         dst.corner[13] = t.corner[17] ;
         dst.corner[15] = t.corner[6] ;
         dst.corner[16] = inc(t.corner[9],1,3) ;
         dst.corner[17] = t.corner[19] ;
         dst.corner[19] = t.corner[10] ;
         dst.center[5] = t.center[8] ;
         dst.center[7] = t.center[9] ;
         dst.center[8] = t.center[11] ;
         dst.center[9] = t.center[5] ;
         dst.center[11] = t.center[7] ;
         break ;
case 1:
         dst.corner[1] = t.corner[13] ;
         dst.corner[2] = t.corner[14] ;
         dst.corner[8] = inc(t.corner[15],1,3) ;
         dst.corner[11] = inc(t.corner[1],1,3) ;
         dst.corner[12] = t.corner[8] ;
         dst.corner[13] = inc(t.corner[2],1,3) ;
         dst.corner[14] = inc(t.corner[11],1,3) ;
         dst.corner[15] = t.corner[17] ;
         dst.corner[17] = t.corner[18] ;
         dst.corner[18] = inc(t.corner[12],2,3) ;
         dst.center[1] = t.center[8] ;
         dst.center[3] = t.center[10] ;
         dst.center[8] = t.center[11] ;
         dst.center[10] = t.center[1] ;
         dst.center[11] = t.center[3] ;
         break ;
case 2:
         dst.corner[2] = t.corner[17] ;
         dst.corner[3] = inc(t.corner[18],1,3) ;
         dst.corner[6] = t.corner[19] ;
         dst.corner[7] = inc(t.corner[11],2,3) ;
         dst.corner[11] = t.corner[13] ;
         dst.corner[13] = inc(t.corner[6],1,3) ;
         dst.corner[14] = inc(t.corner[2],1,3) ;
         dst.corner[17] = t.corner[3] ;
         dst.corner[18] = inc(t.corner[14],1,3) ;
         dst.corner[19] = t.corner[7] ;
         dst.center[2] = t.center[11] ;
         dst.center[4] = t.center[10] ;
         dst.center[7] = t.center[4] ;
         dst.center[10] = t.center[2] ;
         dst.center[11] = t.center[7] ;
         break ;
case 3:
         dst.corner[3] = t.corner[18] ;
         dst.corner[4] = inc(t.corner[19],2,3) ;
         dst.corner[7] = t.corner[4] ;
         dst.corner[10] = inc(t.corner[16],2,3) ;
         dst.corner[12] = t.corner[10] ;
         dst.corner[14] = inc(t.corner[12],1,3) ;
         dst.corner[16] = t.corner[17] ;
         dst.corner[17] = t.corner[14] ;
         dst.corner[18] = inc(t.corner[7],2,3) ;
         dst.corner[19] = inc(t.corner[3],2,3) ;
         dst.center[0] = t.center[9] ;
         dst.center[3] = t.center[10] ;
         dst.center[7] = t.center[3] ;
         dst.center[9] = t.center[7] ;
         dst.center[10] = t.center[0] ;
         break ;
case 4:
         dst.corner[1] = t.corner[5] ;
         dst.corner[2] = t.corner[13] ;
         dst.corner[3] = inc(t.corner[14],1,3) ;
         dst.corner[5] = inc(t.corner[16],1,3) ;
         dst.corner[6] = t.corner[17] ;
         dst.corner[13] = inc(t.corner[15],1,3) ;
         dst.corner[14] = inc(t.corner[1],1,3) ;
         dst.corner[15] = inc(t.corner[6],1,3) ;
         dst.corner[16] = t.corner[3] ;
         dst.corner[17] = inc(t.corner[2],1,3) ;
         dst.center[2] = t.center[8] ;
         dst.center[3] = t.center[2] ;
         dst.center[6] = t.center[7] ;
         dst.center[7] = t.center[3] ;
         dst.center[8] = t.center[6] ;
         break ;
case 5:
         dst.corner[2] = inc(t.corner[18],2,3) ;
         dst.corner[3] = t.corner[19] ;
         dst.corner[4] = t.corner[9] ;
         dst.corner[6] = inc(t.corner[17],2,3) ;
         dst.corner[9] = inc(t.corner[15],1,3) ;
         dst.corner[15] = t.corner[2] ;
         dst.corner[16] = t.corner[6] ;
         dst.corner[17] = inc(t.corner[3],2,3) ;
         dst.corner[18] = t.corner[4] ;
         dst.corner[19] = inc(t.corner[16],2,3) ;
         dst.center[3] = t.center[4] ;
         dst.center[4] = t.center[9] ;
         dst.center[6] = t.center[11] ;
         dst.center[9] = t.center[6] ;
         dst.center[11] = t.center[3] ;
         break ;
      }
      t = dst ;
   }
}
int NMOVES = 6 ;
#include <iostream>
#include <vector>
#include <set>
using namespace std ;
puz solved ;
int maxn = -1 ;
int stack[100] ;
int globald ;
void dfs(const puz &a, int lastm, int togo) {
   if (togo <= 1) {
      int cornr = 0 ;
      int centr = 0 ;
      for (int i=0; i<12; i++)
         if (a.center[i] == solved.center[i])
            centr++ ;
      for (int i=0; i<20; i++)
         if (a.corner[i] == solved.corner[i])
            cornr++ ;
      if (togo == 0 && cornr >= 17 && centr >= 9) {
         cout << " " << cornr << " " << centr << ":" ;
         for (int i=globald; i>0; i--)
            cout << " " << (char)('A'+stack[i]/10) << (stack[i]%10) ;
         cout << endl << flush ;
      }
      if (togo == 0 || (togo == 1 && (cornr < 7 || centr < 3)))
         return ;
   }
   puz b ;
   for (int m=0; m<NMOVES; m++)
      if (m != lastm)
         for (int n=1; n<maxn; n++) {
            stack[togo] = 10 * m + n ;
            slowmove(a, b, m, n) ;
            dfs(b, m, togo-1) ;
         }
}
#include <sys/time.h>
static double start ;
double walltime() {
   struct timeval tv ;
   gettimeofday(&tv, 0) ;
   return tv.tv_sec + 0.000001 * tv.tv_usec ;
}
double duration() {
   double now = walltime() ;
   double r = now - start ;
   start = now ;
   return r ;
}
int main() {
   puz a, b ;
   init(a) ;
   for (int m=0; m<NMOVES; m++) {
      for (int n=1; ; n++) {
         slowmove(a, b, m, n) ;
         if (a == b) {
            if (maxn < 0)
               maxn = n ;
            else if (maxn != n) {
               cout << "Bad cycles; " << n << " vs " << maxn << endl ;
               exit(10) ;
            }
            break ;
         }
      }
   }
   solved = a ;
   duration() ;
   for (int d=0; ; d++) {
      globald = d ;
      dfs(a, -1, d) ;
      cout << "Searched " << d << " in " << duration() << endl << flush ;
   }
}
