#
#   Read a ksolve file (with some significant restrictions) and write out
#   efficient C++ code to simulate that puzzle.
#
#   We only read Set, Solved, and Move blocks for now.  Anything else
#   kills us.
#
my $linebuf = undef ;
sub getline {
   if (defined($linebuf)) {
      my $r = $linebuf ;
      $linebuf = undef ;
      return $r ;
   }
   while (1) {
      my $lin = <> ;
      return undef if !defined($lin) ;
      chomp $lin ;
      $lin =~ s/^\s+// ;
      $lin =~ s/\s+$// ;
      next if $lin =~ /^#/ ;
      next if $lin =~ /^\s*$/ ;
      return $lin ;
   }
}
sub ungetline {
   die "Can't unget multiple lines" if defined($linebuf) ;
   $linebuf = shift ;
}
my @sets = () ;
my %n, %o ;
my @moves = () ;
my $solved ;
my %moves ;
sub readperm {
   my %p = () ;
   for (@sets) {
      $p{$_}[0] = [1..$n{$_}] ;
      $p{$_}[1] = [(0) x $n{$_}] ;
   }
   while (1) {
      my $lin = getline() ;
      last if $lin eq 'End' ;
      $set = $lin ;
      die "Bad perm: undefined [$set]" if !defined($n{$set}) ;
      $lin = getline() ;
      die "Bad permutation [$lin]" if $lin !~ /^[\d\s]+$/ ;
      @pp = map { - - $_ } split " ", $lin ;
      $p{$set}[0] = [@pp] ;
      $lin = getline() ;
      if ($lin =~ /^[\d\s]+$/) {
         @pp = map { - - $_ } split " ", $lin ;
         $p{$set}[1] = [@pp] ;
      } else {
         ungetline($lin) ;
      }
   }
   return \%p ;
}
while (defined($lin = getline())) {
   if ($lin =~ /Set (\S+) (\d+) ?(\d*)/) {
      my $name = $1 ;
      my $n = $2 ;
      my $o = $3 ;
      $o = 1 if !defined($o) || $o eq '' ;
      $n{$name} = $n ;
      $o{$name} = $o ;
      push @sets, $name ;
   } elsif ($lin =~ /Solved/) {
      $solved = readperm() ;
   } elsif ($lin =~ /Move (\S+)/) {
      push @moves, $1 ;
      $moves{$1} = readperm() ;
   } else {
      die "Bad input at [$lin]" ;
   }
}
#
#   For now we write just the permutations and use slow C++ code to
#   actually do the work.  We'll improve this as we go.  We assume the
#   set names are legal C++ identifiers and we use these to define
#   fields of a class structure.  We do make the names be lower case
#   to make things prettier.
#
print <<EOF ;
#include <string.h>
struct puz {
EOF
for (@sets) {
   my $siz = $o{$_} * $n{$_} ;
   die "Puzzle too big; extend range" if $siz > 256 ;
   my $f = lc $_ ;
   print "   char ${f}[$n{$_}] ;\n" ;
}
print <<EOF ;
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
EOF
print "void init(puz &p) {\n" ;
print "   p = puz({\n" ;
for (@sets) {
   my @n = @{$solved->{$_}[0]} ;
   my @o = @{$solved->{$_}[1]} ;
   my $o = $o{$_} ;
   print "         {" ;
   for (my $i=0; $i<$n{$_}; $i++) {
      print "," if $i ;
      print $o * ($n[$i]-1) + $o[$i] ;
   }
   print "},\n" ;
}
print <<EOF ;
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
EOF
my $mvcnt = 0 ;
for $move (@moves) {
   print "case $mvcnt:\n" ;
   $mvcnt++ ;
   for $set (@sets) {
      my $f = lc $set ;
      my @n = @{$moves{$move}{$set}[0]} ;
      my @o = @{$moves{$move}{$set}[1]} ;
      for ($i=0; $i<$n{$set}; $i++) {
         my $o = $o[$n[$i]-1] ;
         next if $o == 0 && $n[$i] == $i+1 ;
         my $ni = $n[$i]-1 ;
         if ($o == 0) {
            print "         dst.${f}[$i] = t.${f}[$ni] ;\n" ;
         } else {
            print "         dst.${f}[$i] = inc(t.${f}[$ni],$o,$o{$set}) ;\n" ;
         }
      }
   }
   print "         break ;\n" ;
}
$nmoves = @moves ;
print <<EOF ;
      }
      t = dst ;
   }
}
int NMOVES = $nmoves ;
#include <iostream>
#include <vector>
#include <set>
using namespace std ;
int main() {
   puz a, b ;
   init(a) ;
   int maxn = -1 ;
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
   vector<puz> q ;
   set<puz> world ;
   q.push_back(a) ;
   world.insert(a) ;
   vector<puz> nextlev ;
   for (int d=0; ; d++) {
      cout << d << " " << q.size() << endl << flush ;
      for (int i=0; i<q.size(); i++)
         for (int m=0; m<NMOVES; m++)
            for (int n=0; n<maxn; n++) {
               slowmove(q[i], b, m, n) ;
               if (world.find(b) == world.end()) {
                  nextlev.push_back(b) ;
                  world.insert(b) ;
               }
            }
      swap(q, nextlev) ;
      nextlev.clear() ;
   }
}
EOF
