#include <iostream>
using namespace std ;
const char *faces = "ABCDEFabcdef" ;
const int NFACES = 12 ;
const int NCORNERS = 20 ;
const char *corners[NCORNERS] = {
   "ABC", "ACD", "ADE", "AEF", "AFB", "BCe", "CDf", "DEb", "EFc", "FBd",
   "bfD", "cbE", "dcF", "edB", "feC", "abc", "acd", "ade", "aef", "afb"
} ;
char symtoface[128] ;
// legal orders.
char illegal[NFACES][NFACES] ;
// how many canonical sequences
const int MAXMOVES = 200 ;
double canonseq[MAXMOVES][NFACES] ;
int main(int argc, char *argv[]) {
   for (int i=0; i<NFACES; i++)
      symtoface[faces[i]] = i ;
   for (int i=0; i<NFACES; i++)
      for (int j=0; j<i; j++)
         illegal[i][j] = 1 ;
   for (int i=0; i<NCORNERS; i++) {
      int jj = 2 ;
      for (int j=0; j<3; j++) {
         illegal[symtoface[corners[i][jj]]][symtoface[corners[i][j]]] = 0 ;
         jj = j ;
      }
   }
   for (int i=0; i<NFACES; i++)
      illegal[i][i] = 1 ;
   cout << 0 << " " << 1 << " " << 1 << endl ;
   for (int i=0; i<NFACES; i++)
      canonseq[1][i] = 4.0 ;
   double total = 1 ;
   for (int mvs=1; mvs+1<MAXMOVES; mvs++) {
      double sum = 0 ;
      for (int i=0; i<NFACES; i++)
         sum += canonseq[mvs][i] ;
      total += sum ;
      cout << mvs << " " << sum << " " << total << endl ;
      for (int i=0; i<NFACES; i++)
         for (int j=0; j<NFACES; j++)
            if (!illegal[i][j])
               canonseq[mvs+1][j] += 4.0 * canonseq[mvs][i] ;
   }
}
