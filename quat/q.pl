#perl
use Math::Trig ;
$|++ ;
#
#   Experiments with quaternions.
#
#  1  1  i  j  k
#  1  1  i  j  k
#  i  i -1  k -j
#  j  j -k -1  i
#  k  k  j -i -1
#
sub mul {
   my $a = shift ;
   my $b = shift ;
   return [$a->[0]*$b->[0]-$a->[1]*$b->[1]-$a->[2]*$b->[2]-$a->[3]*$b->[3],
           $a->[0]*$b->[1]+$a->[1]*$b->[0]+$a->[2]*$b->[3]-$a->[3]*$b->[2],
           $a->[0]*$b->[2]-$a->[1]*$b->[3]+$a->[2]*$b->[0]+$a->[3]*$b->[1],
           $a->[0]*$b->[3]+$a->[1]*$b->[2]-$a->[2]*$b->[1]+$a->[3]*$b->[0]] ;
}
#
#   Now we experiment with what rotations can be generated from a given set
#   of rotations.  Our check of whether our quats work or not is, are the
#   sets so generated of the appropriate size.  In order for this to work,
#   we need to be able to check the distance between quats.
#
sub s2 {
   my $a = shift ;
   my $b = shift ;
   return ($a - $b) * ($a - $b) ;
}
sub d {
   my $a = shift ;
   my $b = shift ;
   return sqrt(s2($a->[0], $b->[0]) + s2($a->[1], $b->[1]) +
               s2($a->[2], $b->[2]) + s2($a->[3], $b->[3])) ;
}
sub cross {
   my $a = shift ;
   my $b = shift ;
   return [0, $a->[2]*$b->[3]-$a->[3]*$b->[2],
           $a->[3]*$b->[1]-$a->[1]*$b->[3], $a->[1]*$b->[2]-$a->[2]*$b->[1]] ;
}
sub dot {
   my $a = shift ;
   my $b = shift ;
   return $a->[1]*$b->[1]+$a->[2]*$b->[2]+$a->[3]*$b->[3] ;
}
sub smul {
   my $a = shift ;
   my $m = shift ;
   return [$a->[0]*$m, $a->[1]*$m, $a->[2]*$m, $a->[3]*$m] ;
}
sub sum {
   my $a = shift ;
   my $b = shift ;
   return [$a->[0]+$b->[0], $a->[1]+$b->[1], $a->[2]+$b->[2], $a->[3]+$b->[3]] ;
}
#
#   We start with cube rotation generators.  I don't see how conjugation
#   by these will retain enough information.
#
my $sp5 = sqrt(0.5) ;
my $s32 = sqrt(3)/2 ;
#
my @g = () ;
#
sub cubical {
   @g = ([$sp5, $sp5, 0, 0], [$sp5, 0, $sp5, 0]) ;
}
sub tetrahedral {
   @g = ([0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, -0.5]) ;
}
sub dodecahedral {
   my $d36 = 2 * pi / 10 ;
   my $dx = 1/2 + 3/10 * sqrt(5) ;
   my $dy = 1/2 + 1/10 * sqrt(5) ;
   my $dd = sqrt($dx * $dx + $dy * $dy) ;
   $dx /= $dd ;
   $dy /= $dd ;
   @g = ([cos $d36, $dx * sin $d36, $dy * sin $d36, 0],
         [0.5, 0.5, 0.5, 0.5]) ;
}
sub icosahedral {
   my $dx = 1/6+sqrt(5)/6 ;
   my $dy = 2/3+sqrt(5)/3 ;
   my $dd = sqrt($dx * $dx + $dy * $dy) ;
   $dx /= $dd ;
   $dy /= $dd ;
   my $ang = 2 * pi / 6 ;
   @g = ([cos $ang, $dx * sin $ang, $dy * sin $ang, 0],
         [cos $ang, -$dx * sin $ang, $dy * sin $ang, 0]) ;
}
sub invrot {
   my $a = shift ;
   return [$a->[0], -$a->[1], -$a->[2], -$a->[3]] ;
}
sub side {
   my $a = shift ;
   return 1 if $a > $eps ;
   return -1 if $a < -$eps ;
   return 0 ;
}
#
#   Q and -Q are the same rotation.  If we are looking at the space we
#   want to choose the one that is positive.  To handle this we use two
#   identities and we divide the state space by two.  This means we do
#   more work, but it saves us from floating point slop.
#
#   To move a position by a rotation, we do a conjugation, but to compose
#   rotations we just do a simple multiplication.  So for now we focus
#   only on the rotations, starting with the null rotation which is
#   [1, 0, 0, 0].
#
my $ident = [1, 0, 0, 0] ;
my $nident = [-1, 0, 0, 0] ;
#
#   Epsilon.
#
$eps = 1e-9 ;
#
#   What symmetry?
#
my $sym = shift ;
if ($sym eq 'cube') {
   cubical() ;
} elsif ($sym eq 'tetrahedron') {
   tetrahedral() ;
} elsif ($sym eq 'dodecahedron') {
   dodecahedral() ;
} elsif ($sym eq 'icosahedron') {
   icosahedral() ;
} else {
   die "Bad shape" ;
}
sub generate {
   my @g = @_ ;
   my @q = ($ident, $nident) ;
   my $i ;
   for (my $i=0; $i<@q; $i++) {
      for (my $j=0; $j<@g; $j++) {
         my $ns = mul($g[$j], $q[$i]) ;
         my $seen = 0 ;
         for (my $k=0; $k<@q; $k++) {
            if (d($ns, $q[$k]) < $eps) {
               $seen++ ;
               last ;
            }
         }
         if (!$seen) {
            push @q, $ns ;
         }
      }
   }
   return @q ;
}
#
#   Solving three planes using linear algebra.  First a 3x3 determinant.
#
sub det3x3 {
   my $a00 = shift ;
   my $a01 = shift ;
   my $a02 = shift ;
   my $a10 = shift ;
   my $a11 = shift ;
   my $a12 = shift ;
   my $a20 = shift ;
   my $a21 = shift ;
   my $a22 = shift ;
   return $a00 * ($a11 * $a22 - $a12 * $a21) +
          $a01 * ($a12 * $a20 - $a10 * $a22) +
          $a02 * ($a10 * $a21 - $a11 * $a20) ;
}
#
#   Find the point at the intersection of three planes, if there is one.
#
sub solvethreeplanes {
   my $pi1 = shift ;
   my $pi2 = shift ;
   my $pi3 = shift ;
   my $planes = shift ;
   my @planes = @{$planes} ;
   my $p1 = $planes[$pi1] ;
   my $p2 = $planes[$pi2] ;
   my $p3 = $planes[$pi3] ;
   my $det = det3x3($p1->[1], $p1->[2], $p1->[3],
                    $p2->[1], $p2->[2], $p2->[3],
                    $p3->[1], $p3->[2], $p3->[3]) ;
   return undef if abs($det) < $eps ;
   my $x = det3x3($p1->[0], $p1->[2], $p1->[3],
                  $p2->[0], $p2->[2], $p2->[3],
                  $p3->[0], $p3->[2], $p3->[3]) / $det ;
   my $y = det3x3($p1->[1], $p1->[0], $p1->[3],
                  $p2->[1], $p2->[0], $p2->[3],
                  $p3->[1], $p3->[0], $p3->[3]) / $det ;
   my $z = det3x3($p1->[1], $p1->[2], $p1->[0],
                  $p2->[1], $p2->[2], $p2->[0],
                  $p3->[1], $p3->[2], $p3->[0]) / $det ;
   my $i ;
   for ($i=0; $i<@planes; $i++) {
      next if $i == $pi1 || $i == $pi2 || $i == $pi3 ;
      my $pt = $planes[$i] ;
      my $dt = $pt->[1] * $x + $pt->[2] * $y + $pt->[3] * $z ;
      return undef if $pt->[0] > 0 && $dt > $pt->[0] ;
      return undef if $pt->[0] < 0 && $dt < $pt->[0] ;
   }
   return [0, $x, $y, $z] ;
}
#
#   Rotate a plane.
#
sub rotateplane {
   my $q = shift ;
   my $p = shift ;
   my $t = mul(mul($q, [0, $p->[1], $p->[2], $p->[3]]), invrot($q)) ;
   $t->[0] = $p->[0] ;
   return $t ;
}
#
#   Find the unique planes.  We actually just generate the subgroup of the
#   rotation group that generates the planes, because from this we can
#   easily regenerate all the planes.
#
sub genuniqueplanes {
   my $q = shift ;
   my $rotations = shift ;
   my @rotations = @{$rotations} ;
   my @planes = () ;
   my @planerot = () ;
   for (my $i=0; $i<@rotations; $i++) {
      my $p = rotateplane($rotations[$i], $q) ;
      my $seen = 0 ;
      for (my $j=0; $j<@planes; $j++) {
         if (d($p, $planes[$j]) <= $eps) {
            $seen++ ;
            last ;
         }
      }
      if (!$seen) {
         push @planes, $p ;
         push @planerot, $rotations[$i] ;
      }
   }
   return @planerot ;
}
#
#   From a set of planes, generate a single face.
#
#   We pick the first plane as the main one.  We then iterate for all other
#   pairs of planes, solve for a vertex, and check that vertex against the
#   other planes.
#
sub getface {
   my @planes = @_ ;
   my @face = () ;
   for (my $i=1; $i<@planes; $i++) {
      for (my $j=$i+1; $j<@planes; $j++) {
         my $p = solvethreeplanes(0, $i, $j, \@planes) ;
         if (defined($p)) {
            my $seen = 0 ;
            for (my $k=0; $k<@face; $k++) {
               if (d($p, $face[$k]) < $eps) {
                  $seen++ ;
                  last ;
               }
            }
            next if $seen ;
            push @face, $p ;
         }
      }
   }
#
#   Sort the points of the face.
#
   while (1) {
      my $changed = 0 ;
      for (my $i=0; $i<@face; $i++) {
         my $j = ($i + 1) % @face ;
         if (dot($planes[0], cross($face[$i], $face[$j])) < 0) {
            my $t = $face[$i] ;
            $face[$i] = $face[$j] ;
            $face[$j] = $t ;
            $changed = 1 ;
         }
      }
      last if !$changed ;
   }
   return @face ;
}
#
#   Normalize:  make the size be 1.
#
sub normalize {
   my $q = shift ;
   my $s = sqrt(dot($q, $q)) ;
   return [$q->[0]/$s, $q->[1]/$s, $q->[2]/$s, $q->[3]/$s] ;
}
#
#   To make a normal from a plane or quat, kill the constant, then
#   normalize.
#
sub makenormal {
   my $q = shift ;
   return normalize([0, $q->[1], $q->[2], $q->[3]]) ;
}
#
#   Get a description of a plane from the command line.
#
sub getplanefromcommandline {
   my $a = shift @ARGV ;
   my $b ;
   my $c ;
   if ($a eq 'face') {
      ($a, $b, $c) = @{$facenormal}[1,2,3] ;
   } elsif ($a eq 'vertex') {
      ($a, $b, $c) = @{$vertexnormal}[1,2,3] ;
   } elsif ($a eq 'edge') {
      ($a, $b, $c) = @{$edgenormal}[1,2,3] ;
   } else {
      $b = shift @ARGV ;
      $c = shift @ARGV ;
   }
   my $d = shift @ARGV ;
   return [$d, $a, $b, $c] ;
}
#
#   First, generate the rotation group.
#
my @rotations = generate(@g) ;
print "// Total is ", scalar @rotations, "\n" ;
#
#   For this rotation group, let's figure out the base set of planes,
#   starting with a plane whose normal is the first generator, and
#   rotating this by the rotation group.  We use this to define the
#   default face normal, vertex normal, and edge normal.
#
my @base = @{$g[0]} ;
my @baseplanerot = genuniqueplanes(\@base, \@rotations) ;
my @baseplanes = map { rotateplane($_, \@base) } @baseplanerot ;
my @baseface = getface(@baseplanes) ;
#
#   From this face we can pick out normals for face, vertex, and edge.
#
$facenormal = makenormal($baseplanes[0]) ;
$edgenormal = makenormal(sum($baseface[0], $baseface[1])) ;
$vertexnormal = makenormal($baseface[0]) ;
print "//Facenormal @{$facenormal}\n" ;
print "//Edgenormal @{$edgenormal}\n" ;
print "//Vertexnormal @{$vertexnormal}\n" ;
#
#   Pull in the actual boundaries.
#
#
my $boundary = getplanefromcommandline() ;
my @planerot = genuniqueplanes($boundary, \@rotations) ;
my @planes = map { rotateplane($_, $boundary) } @planerot ;
$nplanes = @planes ;
print "// Total planes is $nplanes\n" ;
my @face = getface(@planes) ;
#
#   Now do the cuts.  We split the face into multiple faces based on the
#   rotations of the cuts.
#
my $cutplane = getplanefromcommandline() ;
my $d = $cutplane->[0] ;
#
my @faces = [@face] ;
for ($i=0; $i<@rotations; $i++) {
   my $q = rotateplane($rotations[$i], $cutplane) ;
   my @nfaces = () ;
   for ($j=0; $j<@faces; $j++) {
      my @face = @{$faces[$j]} ;
      my @inout = map { side(dot($_, $q) - $d) } @face ;
      my $seen = 0 ;
      for (@inout) {
         $seen |= (1 << ($_ + 1)) ;
      }
      if (($seen & 5) == 5) { # saw both sides
         for (my $s=-1; $s <= 1; $s += 2) {
            my @nface = () ;
            for ($k=0; $k<@face; $k++) {
               if ($inout[$k] == $s || $inout[$k] == 0) {
                  push @nface, $face[$k] ;
               }
               my $kk = ($k + 1) % @face ;
               if ($inout[$k] + $inout[$kk] == 0 && $inout[$k] != 0) {
                  my $vk = dot($face[$k], $q) - $d ;
                  my $vkk = dot($face[$kk], $q) - $d ;
                  my $r = $vk / ($vk - $vkk) ;
                  my $pt = sum(smul($face[$k], (1-$r)), smul($face[$kk], $r)) ;
                  push @nface, $pt ;
               }
            }
            push @nfaces, [@nface] ;
         }
      } else { # no cut
         push @nfaces, [@face] ;
      }
   }
   @faces = @nfaces ;
}
print "// Total faces now is ", scalar @faces, "\n" ;
#
#   Print the result.
#
print "f([\n" ;
for ($i=0; $i<@planerot; $i++) {
   for ($k=0; $k<@faces; $k++) {
      print " [" ;
      my @face = @{$faces[$k]} ;
      for ($j=0; $j<@face; $j++) {
         my $q = rotateplane($planerot[$i], $face[$j]) ; # really point
         print "[$q->[1],$q->[2],$q->[3]]," ;
      }
      print "],\n" ;
   }
}
print "]);\n" ;
