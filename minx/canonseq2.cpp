#include <iostream>
using namespace std ;
const char *faces = "ABCDEFabcdef" ;
const int NFACES = 12 ;
const int NCORNERS = 20 ;
const int NSTATES = (1<<NFACES) ;
const char *corners[NCORNERS] = {
   "ABC", "ACD", "ADE", "AEF", "AFB", "CBe", "DCf", "EDb", "FEc", "BFd",
   "fbD", "bcE", "cdF", "deB", "efC", "acb", "adc", "aed", "afe", "abf"
} ;
char symtoface[128] ;
// legal orders.
int commutes[NFACES] ;
int adjacent[NFACES] ;
// how many canonical sequences
const int MAXMOVES = 200 ;
double canonseq[MAXMOVES][NSTATES] ;
int main(int argc, char *argv[]) {
   for (int i=0; i<NFACES; i++) {
      symtoface[faces[i]] = i ;
      commutes[i] = (1<<NFACES)-1-(1<<i) ;
   }
   for (int i=0; i<NCORNERS; i++) {
      int jj = 2 ;
      for (int j=0; j<3; j++) {
         int f1 = symtoface[corners[i][jj]] ;
         int f2 = symtoface[corners[i][j]] ;
         commutes[f1] &= ~(1<<f2) ;
         jj = j ;
      }
   }
   for (int i=0; i<NFACES; i++)
      commutes[i] &= ~(1<<i) ;
   cout << 0 << " " << 1 << " " << 1 << endl ;
   canonseq[0][0] = 1 ;
   double total = 0 ;
   double osum = 1 ;
   for (int mvs=0; mvs+1<MAXMOVES; mvs++) {
      double sum = 0 ;
      for (int i=0; i<NSTATES; i++) {
         if (canonseq[mvs][i] == 0)
            continue ;
 cout << " " << i << " " << canonseq[mvs][i] << endl ;
         sum += canonseq[mvs][i] ;
      }
      cout.precision(15) ;
      total += sum ;
      cout << mvs << " " << sum << " " << total << " " << (sum / osum) << endl ;
      osum = sum ;
      for (int i=0; i<NSTATES; i++)
         for (int j=0; j<NFACES; j++) {
            if (((i >> j) & 1) == 0 && (i & commutes[j] & ((1<<j)-1)) == 0) {
               int ni = (i & commutes[j]) | (1<<j) ;
               canonseq[mvs+1][ni] += 4.0 * canonseq[mvs][i] ;
            }
         }
   }
}
