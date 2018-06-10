#perl
use Math::Trig ;
$|++ ;
#
#   Experiments with quaternions.  Touch.
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
sub quatangle {
   my $q = shift ;
   return 2*acos($q->[0]) ;
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
#
#   Epsilon.
#
$eps = 1e-9 ;
#
#   What symmetry?
#
my $sym = shift ;
if ($sym eq 'cube' || $sym eq 'c') {
   cubical() ;
} elsif ($sym eq 'tetrahedron' || $sym eq 't') {
   tetrahedral() ;
} elsif ($sym eq 'dodecahedron' || $sym eq 'd') {
   dodecahedral() ;
} elsif ($sym eq 'icosahedron' || $sym eq 'i') {
   icosahedral() ;
} else {
   die "Bad shape" ;
}
for (@g) {
   print ">> [@{$_}]\n" ;
}
sub generate {
   my @g = @_ ;
   my @q = ($ident) ;
   my $i ;
   for (my $i=0; $i<@q; $i++) {
      for (my $j=0; $j<@g; $j++) {
         my $ns = mul($g[$j], $q[$i]) ;
         my $negns = smul($ns, -1) ;
         my $seen = 0 ;
         for (my $k=0; $k<@q; $k++) {
            if (d($ns, $q[$k]) < $eps ||
                d($negns, $q[$k]) < $eps) {
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
#   To normalize a plane, normalize the last three values and fix up d
#   appropriately.  We do not normalize the sign of d.
#
sub normalizeplane {
   my $q = shift ;
   my $s = sqrt($q->[1]*$q->[1]+$q->[2]*$q->[2]+$q->[3]*$q->[3]) ;
   return [$q->[0]/$s, $q->[1]/$s, $q->[2]/$s, $q->[3]/$s] ;
}
#
#   Get a description of a plane from the command line.
#
sub getplanefromcommandline {
   return undef if @ARGV < 2 ;
   my $a = shift @ARGV ;
   my $b ;
   my $c ;
   if ($a eq 'face' || $a eq 'f') {
      ($a, $b, $c) = @{$facenormal}[1,2,3] ;
   } elsif ($a eq 'vertex' || $a eq 'v') {
      ($a, $b, $c) = @{$vertexnormal}[1,2,3] ;
   } elsif ($a eq 'edge' || $a eq 'e') {
      ($a, $b, $c) = @{$edgenormal}[1,2,3] ;
   } else {
      if (@ARGV < 3) {
         unshift @ARGV, $a ;
         return undef ;
      }
      $b = shift @ARGV ;
      $c = shift @ARGV ;
   }
   my $d = shift @ARGV ;
   $d = eval($d) ;
   return [$d, $a, $b, $c] ;
}
#
#   Given a set of faces, cut them by a plane and return a new set of faces.
#
sub cutfaces {
   my $q = shift ;
   my $d = $q->[0] ;
   my @faces = @_ ;
   my @nfaces = () ;
   for (my $j=0; $j<@faces; $j++) {
      my @face = @{$faces[$j]} ;
      my @inout = map { side(dot($_, $q) - $d) } @face ;
      my $seen = 0 ;
      for (@inout) {
         $seen |= (1 << ($_ + 1)) ;
      }
      if (($seen & 5) == 5) { # saw both sides
         for (my $s=-1; $s <= 1; $s += 2) {
            my @nface = () ;
            for (my $k=0; $k<@face; $k++) {
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
   return @nfaces ;
}
#
#   Which side of a plane is a face on?  Assume it is not cut by the plane.
#
sub faceside {
   my $plane = shift ;
   my @face = @_ ;
   my $d = $plane->[0] ;
   for (@face) {
      my $s = side(dot($_, $plane) - $d) ;
      return $s if $s ;
   }
   die "Could not place the face on a side of the plane" ;
}
#
#   Expand the faces by the rotation matrix.
#
sub expandfaces {
   my $planerot = shift ;
   my $faces = shift ;
   my @planerot = @{$planerot} ;
   my @faces = @{$faces} ;
   my @nfaces = () ;
   for (my $i=0; $i<@planerot; $i++) {
      for (my $k=0; $k<@faces; $k++) {
         my @face = @{$faces[$k]} ;
         my @nface = () ;
         for (my $j=0; $j<@face; $j++) {
            my $q = rotateplane($planerot[$i], $face[$j]) ; # really point
            push @nface, $q ;
         }
         push @nfaces, [@nface] ;
      }
   }
   return @nfaces ;
}
#
#   Are two planes the same?  Take into account the fact that the
#   move normals might be negated, or that they may be multiples of
#   each other.
#
sub sameplane {
   my $a = normalize(shift) ;
   my $b = normalize(shift) ;
   return 1 if d($a, $b) < $eps ;
   return 1 if d($a, smul($b, -1)) < $eps ;
   return 0 ;
}
#
#   Find a center of mass for a face, which we use as a key.
#
sub centermassface {
   my $face = shift ;
   my $s = [0, 0, 0, 0] ;
   my $n = 0 ;
   for my $pt (@{$face}) {
      $s = sum($pt, $s) ;
      $n++ ;
   }
   return smul($s, 1/$n) ;
}
#
#   Find a center of mass for a cubie, which we use as a key.
#
sub centermasscubie {
   my $cubie = shift ;
   my $s = [0, 0, 0, 0] ;
   my $n = 0 ;
   for my $face (@{$cubie}) {
      for my $pt (@{$face}) {
         $s = sum($pt, $s) ;
         $n++ ;
      }
   }
   return smul($s, 1/$n) ;
}
#
#   We can find a cubie from the center of mass by the point sets.
#
my @moveplanesets = () ;
#
#   Find a cubie by checking its plane.  We find how many planes in each
#   plane set the cubie is "in".
#
sub keyface {
   my $face = shift ;
   my @s = "" ;
   for my $mplaneset (@moveplanesets) {
      my $t = 0 ;
      for my $mplane (@{$mplaneset}) {
         $t++ if faceside($mplane, @{$face}) > 0 ;
      }
      push @s, $t ;
   }
   return "@s" ;
}
#
#   Where we look up cubie indicies.
#
my %cubiekey ;
my %facelists ;
#
#   Find a cubie, from the key face.
#
sub findcubie {
   my $cubie = shift ;
   my $key = keyface($cubie->[0]) ;
   die "Miss?" if !defined($cubiekey{$key}) ;
   return $cubiekey{$key} ;
}
#
#   Same face (by center of mass)
#
sub sameface {
   my $f1 = shift ;
   my $f2 = shift ;
   return abs(d(centermassface($f1), centermassface($f2))) < $eps ;
}
#
#   Find a face.
#
sub findface {
   my $face = shift ;
   my $cm = centermassface($face) ;
   my $key = keyface($face) ;
   for my $face2 (@{$facelists{$key}}) {
      if (abs(d($cm, centermassface($faces[$face2]))) < $eps) {
         return $face2 ;
      }
   }
   die "Could not find face" ;
}
#
#   Rotate things.
#
sub rotatepoint {
   my $q = shift ;
   my $p = shift ;
   return mul(mul($q, $p), invrot($q)) ;
}
#
sub rotateface {
   my $q = shift ;
   return map { rotatepoint($q, $_) } @_ ;
}
#
sub rotatecubie {
   my $q = shift ;
   my $cubie = shift ;
   return [map { [rotateface($q, @{$_})] } @{$cubie}] ;
}
#
#   Print all the faces.
#
sub showf {
   my @faces = @_ ;
   print "f([\n" ;
   for (my $k=0; $k<@faces; $k++) {
      print " [" ;
      my @face = @{$faces[$k]} ;
      for (my $j=0; $j<@face; $j++) {
         my $q = $face[$j] ;
         print "[$q->[1],$q->[2],$q->[3]]," ;
      }
      print "],\n" ;
   }
   print "]);\n" ;
}
#
#   Comment character.
#
my $comment = '//' ;
for (@ARGV) {
   $comment = '#' if $_ eq 'ksolve' || $_ eq 'gap' ;
}
print "$comment @ARGV\n" ;
#
#   First, generate the rotation group.
#
my @rotations = generate(@g) ;
print "$comment Total is ", scalar @rotations, "\n" ;
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
print "Base face has ", scalar @baseface, " vertices.\n" ;
#
#   From this face we can pick out normals for face, vertex, and edge.
#
$facenormal = makenormal($baseplanes[0]) ;
print "Baseface 0 [@{$baseface[0]}]\n" ;
print "Baseface 1 [@{$baseface[1]}]\n" ;
$edgenormal = makenormal(sum($baseface[0], $baseface[1])) ;
$vertexnormal = makenormal($baseface[0]) ;
print "$comment Facenormal @{$facenormal}\n" ;
print "$comment Edgenormal @{$edgenormal}\n" ;
print "$comment Vertexnormal @{$vertexnormal}\n" ;
#
#   Pull in the actual boundaries.
#
my $boundary = [@{$facenormal}] ;
$boundary->[0] = 1 ;
if ($ARGV[0] eq 'boundary') {
   $boundary = getplanefromcommandline() ;
}
my @planerot = genuniqueplanes($boundary, \@rotations) ;
my @planes = map { rotateplane($_, $boundary) } @planerot ;
$nplanes = @planes ;
print "$comment Total planes is $nplanes\n" ;
my @face = getface(@planes) ;
#
#   Let's show the face here.
#
for ($i=0; $i<@face; $i++) {
#  print "# PT $face[$i][0] $face[$i][1] $face[$i][2] $face[$i][3]\n" ;
}
#
#   Now do the cuts.  We split the face into multiple faces based on the
#   rotations of the cuts.
#
my @cutplanes = () ;
while (@ARGV) {
   my $cutplane = getplanefromcommandline() ;
   last if !defined($cutplane) ;
   push @cutplanes, $cutplane ;
}
#
@faces = [@face] ;
my @moveplanes = () ;
for my $cutplane (@cutplanes) {
   for ($i=0; $i<@rotations; $i++) {
      my $q = rotateplane($rotations[$i], $cutplane) ;
      $seen = 0 ;
      for ($j=0; $j<@moveplanes; $j++) {
         if (sameplane($q, $moveplanes[$j])) {
            $seen++ ;
            last ;
         }
      }
      if ($seen == 0) {
         push @moveplanes, $q ;
         @faces = cutfaces($q, @faces) ;
      }
   }
}
print "$comment Faces now is ", scalar @faces, "\n" ;
print "$comment Unique move planes is ", scalar @moveplanes, "\n" ;
@faces = expandfaces(\@planerot, \@faces) ;
print "$comment Final total faces now is ", scalar @faces, "\n" ;
#
#   Split moveplanes into a list of parallel planes.
#
for ($i=0; $i<@moveplanes; $i++) {
   my $seen = 0 ;
   my $q = $moveplanes[$i] ;
   my $qnormal = makenormal($q) ;
   for ($j=0; $j<@moveplanesets; $j++) {
      if (sameplane($qnormal, makenormal($moveplanesets[$j][0]))) {
         push @{$moveplanesets[$j]}, $q ;
         $seen++ ;
         last ;
      }
   }
   if ($seen == 0) {
      push @moveplanesets, [$q] ;
   }
}
#
#   Make all the parallel planes face the same way, and then sort by d.
#
for my $moveplaneset (@moveplanesets) {
   my @a = map { normalizeplane($_) } @{$moveplaneset} ;
   my $goodnormal = makenormal($a[0]) ;
   for (my $i=0; $i<@a; $i++) {
      if (d(makenormal($a[$i]), $goodnormal) > $eps) {
         $a[$i] = [-$a[$i][0], -$a[$i][1], -$a[$i][2], -$a[$i][3]] ;
      }
   }
   @a = sort { $a->[0] <=> $b->[0] } @a ;
   $moveplaneset = [@a] ;
}
my @sizes = map { scalar @{$_} } @moveplanesets ;
print "$comment move plane sets: [@sizes]\n" ;
#
#   For each set of move planes, find the rotations that are relevant.
#
my @moverotations = () ;
for ($i=0; $i<@rotations; $i++) {
   my $seen = 0 ;
   my $q = $rotations[$i] ;
   next if abs(abs($q->[0])-1) < $eps ;
   my $qnormal = makenormal($q) ;
   for ($j=0; $j<@moveplanesets; $j++) {
      if (sameplane($qnormal, makenormal($moveplanesets[$j][0]))) {
         push @{$moverotations[$j]}, $q ;
         $seen++ ;
         last ;
      }
   }
}
#
#   Sort the rotations by the angle of the rotation.  A bit tricky because
#   while the norms should be the same, they need not be.  So we start by
#   making the norms the same.
#
for my $moverotationlist (@moverotations) {
   my @a = @{$moverotationlist} ;
   my $goodnormal = makenormal($a[0]) ;
   for (my $i=0; $i<@a; $i++) {
      if (d(makenormal($a[$i]), $goodnormal) > $eps) {
         $a[$i] = [-$a[$i][0], -$a[$i][1], -$a[$i][2], -$a[$i][3]] ;
      }
   }
   my @angles = sort { $a->[0] <=> $b->[0] } map { [quatangle($_), $_] } @a ;
   $moverotationlist = [map { $_->[1] } @angles] ;
}
@sizes = map { scalar @{$_} } @moverotations ;
print "$comment move rotation sets: [@sizes]\n" ;
#
#   Now, break the faces up into cubie sets according to which side of
#   each plane each is on.  Each should be on a single side of a
#   cutting plane (and not entirely on any particular plane).  We build
#   a string describing sides and collect faces.
#
my %cubies ;
for (my $i=0; $i<@faces; $i++) {
   my $face = $faces[$i] ;
   my $s = keyface($face) ;
   push @{$facelists{$s}}, $i ;
   push @{$cubies{$s}}, $face ;
}
print "$comment Count of cubies is ", scalar keys %cubies, "\n" ;
#
#   Each cubie has a "center of mass" which we approximate as the
#   sum of all the points in all the faces (counting ones multiply that
#   appear multiple times).  We use this center of mass to "identify"
#   a cubie, and with this we will calculate the orbit of all the
#   cubies.  For now we do not "hold anything still"; we will deal with
#   that later.
#
my @cubies = () ;
my @cubiekeys = () ;
for my $key (keys %cubies) {
   $cubiekey{$key} = scalar @cubies ;
   push @cubies, $cubies{$key} ;
   push @cubiekeys, $key ;
}
#
#   Sort the cubies around each *corner* so they are clockwise.  Only
#   relevant for cubies that have more than two faces.  Note that there
#   may be many points; vertex corners on an icosahedron have five
#   faces, and of course in general for convex shapes there is no bound
#   on the maximum number of faces in a corner.
#
for (my $i=0; $i<@cubies; $i++) {
   my @cubie = @{$cubies[$i]} ;
   next if @cubie < 3 ;
   my $s = keyface($cubie[0]) ;
   my @facelist = @{$facelists{$s}} ;
   my @cm = map { centermassface($_) } @cubie ;
   my $cm = centermassface([@cm]) ;
   while (1) {
      my $changed = 0 ;
      for (my $i=0; $i<@cubie; $i++) {
         my $j = ($i + 1) % @cubie ;
         my $v = dot($cm, cross($cm[$i], $cm[$j])) ;
         if (dot($cm, cross($cm[$i], $cm[$j])) < 0) {
            my $t = $cubie[$i] ;
            $cubie[$i] = $cubie[$j] ;
            $cubie[$j] = $t ;
            $t = $cm[$i] ;
            $cm[$i] = $cm[$j] ;
            $cm[$j] = $t ;
            $t = $facelist[$i] ;
            $facelist[$i] = $facelist[$j] ;
            $facelist[$j] = $t ;
            $changed = 1 ;
         }
      }
      last if !$changed ;
   }
   $cubies[$i] = [@cubie] ;
   $facelists{$s} = [@facelist] ;
}
#
for my $moverotationlist (@moverotations) {
   my @a = @{$moverotationlist} ;
   my $goodnormal = makenormal($a[0]) ;
   for (my $i=0; $i<@a; $i++) {
      if (d(makenormal($a[$i]), $goodnormal) > $eps) {
         $a[$i] = [-$a[$i][0], -$a[$i][1], -$a[$i][2], -$a[$i][3]] ;
      }
   }
   my @angles = sort { $a->[0] <=> $b->[0] } map { [quatangle($_), $_] } @a ;
   $moverotationlist = [map { $_->[1] } @angles] ;
}
#
#   Build an array that takes each face to a cubie ordinal and a face
#   number.
#
my @facetocubies ;
for (my $i=0; $i<@faces; $i++) {
   my $key = keyface($faces[$i]) ;
   for (my $j=0; $j<@{$facelists{$key}}; $j++) {
      if ($i == $facelists{$key}[$j]) {
         push @facetocubies, [$cubiekey{$key}, $j] ;
         last ;
      }
   }
}
#
#   Now we do a breadth-first search from each unseen cubie calculating the
#   orbits.
#
my @typename = qw(? CENTER EDGE CORNER C4RNER C5RNER C6RNER C7RNER ? ? ? ? ?) ;
my @cubiesetname = () ;
my @cubietypecounts = () ;
my @orbitoris = () ;
my @seen = () ;
my $cubiesetnum = -1 ;
my @cubiesetnum ;
my @cubieordnum ;
my @cubieord = () ;
for (my $i=0; $i<@cubies; $i++) {
   next if $seen[$i] ;
   my $cubie = $cubies[$i] ;
   $cubiesetnum++ ;
   $cubiesetnum[$i] = $cubiesetnum ;
   my $facecnt = scalar @{$cubie} ;
   my $typectr = 0+$cubietypecounts[$facecnt]++ ;
   my $typename = $typename[$facecnt] . ($typectr == 0 ? '' : ($typectr+1)) ;
   $cubiesetname[$cubiesetnum] = $typename ;
   $orbitoris[$cubiesetnum] = $facecnt ;
   my @q = ($i) ;
   my $qg = 0 ;
   $seen[$i]++ ;
   while ($qg < @q) {
      $s = $q[$qg++] ;
      $cubiesetnum[$s] = $cubiesetnum ;
      $cubieordnum[$s] = 0+$cubieord[$cubiesetnum]++ ;
      for my $movesets (@moverotations) {
         for my $rotation (@{$movesets}) {
            my $tq = findcubie(rotatecubie($rotation, $cubies[$s])) ;
            next if $seen[$tq] ;
            push @q, $tq ;
            $seen[$tq]++ ;
         }
      }
   }
}
print "$comment Cubie sets are [@cubiesetnum]\n" ;
#
#   Test face twists.
#
my $mvcnt = 0 ;
my @movesbyslice = () ;
my @cmovesbyslice = () ;
for (my $k=0; $k<@moveplanesets; $k++) {
   my @moveplaneset = @{$moveplanesets[$k]} ;
   my @slicenum = () ;
   my @slicecnts = () ;
   for (my $i=0; $i<@faces; $i++) {
      my $face = $faces[$i] ;
      my $t = 0 ;
      for my $moveplane (@moveplaneset) {
         $t++ if faceside($moveplane, @{$face}) > 0 ;
      }
      $slicenum[$i] = $t ;
      $slicecnts[$t]++ ;
   }
   print "$comment Slicecounts are [@slicecnts]\n" ;
   # do moves; single slice moves.
   my @axismoves = () ;
   my @axiscmoves = () ;
   for (my $sc=0; $sc<@slicecnts; $sc++) {
      my $mv = '' ;
      my @slicemoves = () ;
      my @slicecmoves = () ;
      my @cubiedone = () ;
      for (my $i=0; $i<@faces; $i++) {
         next if $slicenum[$i] != $sc ;
         my @a = ($i) ;
         my $cubie = $facetocubies[$i][0] ;
         my $ori = $facetocubies[$i][1] ;
         my @b = @{$facetocubies[$i]} ;
         my $face = $faces[$i] ;
         my $fi2 = $i ;
         while (1) {
            $slicenum[$fi2] = -1 ;
            my $face2 = [rotateface($moverotations[$k][0], @{$face})] ;
            $fi2 = findface($face2) ;
            last if $slicenum[$fi2] != $sc ;
            push @a, $fi2 ;
            push @b, @{$facetocubies[$fi2]} ;
            $face = $face2 ;
         }
         push @slicemoves, [@a] if @a > 1 ;
         push @slicecmoves, [@b] if @b > 2 && !$cubiedone[$b[0]] ;
         for (my $j=0; $j<@b; $j += 2) {
            $cubiedone[$b[$j]]++ 
         }
      }
      push @axismoves, [@slicemoves] ;
      push @axiscmoves, [@slicecmoves] ;
   }
   push @movesbyslice, [@axismoves] ;
   push @cmovesbyslice, [@axiscmoves] ;
}
#
#   The sets of bitmasks giving moves for each axis based on the number
#   of slices in each axis.  For now only OBTM.
#
my @movesets = ([], [],
   [1], [1, 4], [1, 3, 7], [1, 3, 16, 24], [1, 3, 7, 15, 31],
   [1, 3, 7, 64, 96, 112]) ;
my $allmoves = 0 ;
#
sub getmovesets {
   my $slices = shift ;
   if ($allmoves) {
      return (1, 2) if $slices == 2 ;
      die "All moves not defined if slices not equal 2" ;
   }
   return @{$movesets[$slices]} ;
}
#
#   Write out a ksolve definition file.
#
sub writeksolve {
   # first figure out what sets move.
   my @setmoves ;
   for (my $k=0; $k<@moveplanesets; $k++) {
      my @moveplaneset = @{$moveplanesets[$k]} ;
      my $slices = 1+scalar @moveplaneset ;
      my @moveset = getmovesets($slices) ;
      my $allbits = 0 ;
      for (my $i=0; $i<@moveset; $i++) {
         $allbits |= $moveset[$i] ;
      }
      die "Bad moveset" if @moveset == 0 ;
      my @axiscmoves = @{$cmovesbyslice[$k]} ;
      for (my $i=0; $i<@axiscmoves; $i++) {
         next if (($allbits >> $i) & 1) == 0 ;
         my @slicecmoves = @{$axiscmoves[$i]} ;
         for (my $j=0; $j<@slicecmoves; $j++) {
            $setmoves[$cubiesetnum[$slicecmoves[$j][0]]]++ ;
         }
      }
   }
   for (my $i=0; $i<@cubiesetname; $i++) {
      next if !$setmoves[$i] ;
      print "Set $cubiesetname[$i] $cubieord[$i] $orbitoris[$i]\n" ;
   }
   print "\n" ;
   print "Solved\n" ;
   for (my $i=0; $i<@cubiesetname; $i++) {
      next if !$setmoves[$i] ;
      print "$cubiesetname[$i]\n" ;
      print join " ", 1..$cubieord[$i], "\n" ;
   }
   print "End\n" ;
   print "\n" ;
   # Iterate through the axis, pick up the move sets we care about,
   # iterate through those.  For now give arbitrary names.
   my $movename = 'A' ;
   for (my $k=0; $k<@moveplanesets; $k++) {
      my @moveplaneset = @{$moveplanesets[$k]} ;
      my $slices = 1+scalar @moveplaneset ;
      my @moveset = getmovesets($slices) ;
      for (my $i=0; $i<@moveset; $i++) {
         my $movebits = $moveset[$i] ;
         print "Move $movename\n" ;
         my @perm = () ;
         my @ori = () ;
         for (my $i=0; $i<@cubiesetname; $i++) {
            push @perm, [0..$cubieord[$i]-1] ;
            push @ori, [(0) x $cubieord[$i]] ;
         }
         $movename++ ;
         my @axiscmoves = @{$cmovesbyslice[$k]} ;
         for (my $i=0; $i<@axiscmoves; $i++) {
            next if (($movebits >> $i) & 1) == 0 ;
            my @slicecmoves = @{$axiscmoves[$i]} ;
            for (my $j=0; $j<@slicecmoves; $j++) {
               my @mperm = @{$slicecmoves[$j]} ;
               my $setnum = $cubiesetnum[$mperm[0]] ;
               for (my $ii=0; $ii<@mperm; $ii += 2) {
                  $mperm[$ii] = $cubieordnum[$mperm[$ii]] ;
               }
# print ">>> [@mperm]\n" ;
               for (my $ii=0; $ii<@mperm; $ii += 2) {
                  $perm[$setnum][$mperm[$ii]] = $mperm[($ii+2)%@mperm] ;
# print " $mperm[($ii+3)%@mperm] from $mperm[($ii+1)%@mperm]\n" ;
                  $ori[$setnum][$mperm[$ii]] =
      ($mperm[($ii+@mperm-1)%@mperm] - $mperm[($ii+1)%@mperm] + $orbitoris[$setnum]) %
                                                      $orbitoris[$setnum] ;
               }
            }
         }
         for (my $i=0; $i<@cubiesetname; $i++) {
            next if !$setmoves[$i] ;
            print "$cubiesetname[$i]\n" ;
            print join(" ", map { 1+$_ } @{$perm[$i]}), "\n" ;
            if ($orbitoris[$i] > 1) {
               print join(" ", @{$ori[$i]}), "\n" ;
            }
         }
         print "End\n" ;
         print "\n" ;
      }
   }
}
#
#   Write out gap permutations.
#
sub writegap {
   my @perms = () ;
   my $movename = "MA" ;
   for (my $k=0; $k<@moveplanesets; $k++) {
      my @moveplaneset = @{$moveplanesets[$k]} ;
      my $slices = 1+scalar @moveplaneset ;
      my @moveset = getmovesets($slices) ;
      for (my $i=0; $i<@moveset; $i++) {
         $order = undef ;
         my $move = "" ;
         my $movebits = $moveset[$i] ;
         my @axismoves = @{$movesbyslice[$k]} ;
         for (my $i=0; $i<@axismoves; $i++) {
            next if (($movebits >> $i) & 1) == 0 ;
            my @slicemoves = @{$axismoves[$i]} ;
            for (my $j=0; $j<@slicemoves; $j++) {
               # gap perms start with 1
               my @mperm = map { 1+$_ } @{$slicemoves[$j]} ;
               $order = scalar @mperm ;
               $move .= "(" . (join ",", @mperm) . ")" ;
            }
         }
         print "$movename:=$move;\n" ;
         push @perms, $movename ;
         for (my $j=2; $j<$order; $j++) {
            push @perms, "$movename^$j" ;
         }
         $movename++ ;
      }
   }
   print "Gen:=[" ;
   print join ",", @perms ;
   print "];\n" ;
}
#
while (@ARGV) {
   my $a = shift @ARGV ;
   if ($a eq 'ksolve') {
      writeksolve() ;
   } elsif ($a eq 'gap') {
      writegap() ;
   } elsif ($a eq 'threejs') {
      showf(@faces) ;
   } elsif ($a eq 'allmoves') {
      $allmoves++ ;
   } else {
      die "What would you like me to do?" ;
   }
}
