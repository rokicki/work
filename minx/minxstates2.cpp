#include <iostream>
#include <vector>
#include <map>
#include <set>
using namespace std ;
// this is 12 x 5 x 5 sets of move perms.
const int NFACES = 12 ;
const int NMOVES = NFACES * 4 ;
unsigned char perms[] = {
   67,78,90,101,112,69,80,93,104,115,75,86,91,102,113,122,123,125,127,129,124,126,128,131,130,
   1,8,6,4,2,3,9,10,7,5,17,28,39,50,59,18,29,40,51,60,21,32,43,54,62,
   34,46,114,123,89,36,49,117,126,91,42,47,115,124,97,100,101,103,105,107,102,104,106,109,108,
   8,61,72,81,26,9,65,76,84,29,10,62,73,82,27,12,19,17,15,13,14,20,21,18,16,
   45,56,68,125,100,47,58,71,128,102,53,64,69,126,108,111,112,114,116,118,113,115,117,120,119,
   6,19,83,94,37,7,21,84,95,38,10,20,87,98,40,23,30,28,26,24,25,31,32,29,27,
   13,63,111,127,79,14,65,119,128,80,16,64,113,131,82,67,74,72,70,68,69,75,76,73,71,
   4,48,105,96,30,5,49,106,98,32,7,51,109,97,31,34,35,37,39,41,36,38,40,43,42,
   12,74,129,92,24,14,75,130,95,27,20,76,131,93,25,78,85,83,81,79,80,86,87,84,82,
   2,57,116,107,41,3,58,117,109,43,5,60,120,108,42,45,46,48,50,52,47,49,51,54,53,
   23,35,103,122,85,25,38,106,124,86,31,36,104,130,87,89,90,92,94,96,91,93,95,98,97,
   1,52,118,70,15,3,53,119,73,18,9,54,120,71,16,56,63,61,59,57,58,64,65,62,60,
} ;
typedef vector<unsigned char> state ;
typedef vector<unsigned char> mve ;
vector<mve> moves ;
state init ;
set<state> world ;
state makemove(const state &pos, const mve &m) {
   state newpos = pos ;
   for (int i=0; i<5; i++) {
      int oj = 4 ;
      for (int j=0; j<5; j++) {
         newpos[m[5*i+oj]] = pos[m[5*i+j]] ;
         oj = j ;
      }
   }
   return newpos ;
}
typedef vector<unsigned char> cstate ;
cstate compress(const state &pos) {
   cstate r(132/3) ;
   for (int i=0; i<44; i++)
      r[i] = pos[3*i]+6*(pos[3*i+1]+6*pos[3*i+2]) ;
   return r ;
}
void recur(const state &pos, int togo) {
   if (togo == 0) {
      cstate cs = compress(pos) ;
      if (world.find(cs) == world.end())
         world.insert(cs) ;
      return ;
   }
   for (int m=0; m<moves.size(); m++)
      recur(makemove(pos, moves[m]), togo-1) ;
}
int main(int argc, char *argv[]) {
   for (int i=0; i<NFACES; i++)
      for (int j=1; j<5; j++) { // twist
         mve newmove ;
         for (int k=0; k<5; k++)
            for (int ii=0; ii<5; ii++) {
               int base = i * 25 + k * 5 ;
               int src = base + (ii * j) % 5 ;
               newmove.push_back(perms[src]) ;
            }
         moves.push_back(newmove) ;
      }
   for (int i=0; i<NFACES; i++)
      for (int j=0; j<11; j++)
         init.push_back(i) ;
   int ostates = 0 ;
   for (int d=0; ; d++) {
      recur(init, d) ;
      cout << "At " << d << " " << (world.size()-ostates) << endl << flush ;
      ostates = world.size() ;
   }
}
