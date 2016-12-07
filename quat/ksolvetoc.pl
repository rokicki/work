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
sub mulone {
   my $set = shift ;
   my $a = shift ;
   my $b = shift ;
   my @c = () ;
   my $n = $n{$set} ;
   my $o = $o{$set} ;
   my $i ;
   for ($i=0; $i<$n{$set}; $i++) {
      $c[0][$i] = $b->[0][$a->[0][$i]-1] ;
      $c[1][$c[0][$i]-1] = ($b->[1][$b->[0][$a->[0][$i]-1]-1] +
                   $a->[1][$a->[0][$i]-1]) % $o{$set} ;
   }
   return [@c] ;
}
sub mulperm {
   my $a = shift ;
   my $b = shift ;
   my %c = {} ;
   for (@sets) {
      $c{$_} = mulone($_, $a->{$_}, $b->{$_}) ;
   }
   return \%c ;
}
sub isone {
   my $set = shift ;
   my $p = shift ;
   my $i ;
   for ($i=0; $i<$n{$set}; $i++) {
      return 0 if $p->[1][$i] ;
      return 0 if $p->[0][$i] != $i+1 ;
   }
   return 1 ;
}
sub isoneperm {
   my $p = shift ;
   for (@sets) {
      return 0 if !isone($_, $p->{$_}) ;
   }
   return 1 ;
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
sub enc {
   my $n = shift ;
   my $o = shift ;
   my $nn = shift ;
   my $oo = shift ;
   return $n if $oo <= 1 ;
   return 2 * $n + $o if $oo == 2 ;
   return $n + $o * $nn ;
}
sub inc {
   my $n = shift ;
   my $o = shift ;
   my $nn = shift ;
   my $oo = shift ;
   return "$n" if $oo <= 1 || $o == 0 ;
   return "$n^1" if $oo == 2 ;
   my $nnoo = $nn * $oo ;
   my $nno = $nn * $o ;
   return "mod$nnoo($n+$nno)" ;
}
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
      print enc($n[$i]-1, $o[$i], $n, $o) ;
   }
   print "},\n" ;
}
print <<EOF ;
       }) ;
}
EOF
my %modseen ;
for $set (@sets) {
   next if $o{$set} < 3 ;
   my $nnoo = $o{$set} * $n{$set} ;
   next if $modseen{$nnoo}++ ;
   print "#define mod${nnoo}(a) mod${nnoo}a[a]\n" ;
   my $mod2 = $nnoo * 2 ;
   print "char mod${nnoo}a[$mod2] ;\n" ;
}
print <<EOF ;
void slowmove(const puz &src, puz &dst, int m) {
   dst = src ;
   switch (m) {
EOF
my $maxn = 0 ;
my $mvcnt = 0 ;
for $move (@moves) {
   my $p = $moves{$move} ;
   my $n = 1 ;
   while (!isoneperm($p)) {
      print "case $mvcnt:\n" ;
      $mvcnt++ ;
      for $set (@sets) {
         my $f = lc $set ;
         my @n = @{$p->{$set}[0]} ;
         my @o = @{$p->{$set}[1]} ;
         for ($i=0; $i<$n{$set}; $i++) {
            my $o = $o[$n[$i]-1] ;
            next if $o == 0 && $n[$i] == $i+1 ;
            my $ni = $n[$i]-1 ;
            print "      dst.${f}[$i] = ", inc("src.${f}[$ni]", $o, $n{$set}, $o{$set}), " ;\n" ;
         }
      }
      print "      break ;\n" ;
      $n++ ;
      $maxn = $n-1 ;
      $p = mulperm($p, $moves{$move}) ;
   }
}
$nmoves = @moves ;
print <<EOF ;
   }
}
const int NMOVES = $nmoves ;
const int NTWIST = $maxn ;
void slowmove(const puz &src, puz &dst, int m, int n) {
   slowmove(src, dst, m*NTWIST+n-1) ;
}
#include <iostream>
#include <vector>
#include <set>
using namespace std ;
puz solved ;
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
         for (int n=1; n<=NTWIST; n++) {
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
EOF
for $mod (keys %modseen) {
   my $mod2 = $mod * 2 ;
   print <<EOF ;
   for (int i=0; i<$mod2; i++)
      mod${mod}a[i] = i % $mod ;
EOF
}
print <<EOF ;
   init(a) ;
   solved = a ;
   duration() ;
   for (int d=0; ; d++) {
      globald = d ;
      dfs(a, -1, d) ;
      cout << "Searched " << d << " in " << duration() << endl << flush ;
   }
}
EOF
