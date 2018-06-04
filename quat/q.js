"use strict" ;

//  Global epsilon; any difference less than this is ignored.
//  We need to package this better.
//
//  This code is *not* efficient in its puzzle construction, as it
//  always uses simple lists and simple distance functions when
//  uniquifying sets of objects.  We accept this to keep the code
//  simple for now.

var eps = 1e-9 ;


// We need a quaternion class.  We use this to represent rotations,
// planes, and points.

function Quat(a_, b_, c_, d_) {
   if (this instanceof Quat) {
      this.a = a_ ; this.b = b_ ; this.c = c_ ; this.d = d_ ;
   } else {
      return new Quat(a_, b_, c_, d_) ;
   }
}
Quat.prototype = {
   'a':0, 'b':0, 'c':0, 'd':0, // quaternion fields; 1 i j k
   'mul': // quaternion multiplication
   function(q) {
      return Quat(
           this.a*q.a-this.b*q.b-this.c*q.c-this.d*q.d,
           this.a*q.b+this.b*q.a+this.c*q.d-this.d*q.c,
           this.a*q.c-this.b*q.d+this.c*q.a+this.d*q.b,
           this.a*q.d+this.b*q.c-this.c*q.b+this.d*q.a) ;
   },
   'toString': // pretty-print a quat
   function() {
      return '[' + this.a + ',' + this.b + ',' + this.c + ',' + this.d + ']' ;
   },
   'dist': // Euclidean distance between two quaternions
   function(q) {
      return Math.hypot(this.a-q.a, this.b-q.b, this.c-q.c, this.d-q.d) ;
   },
   'cross': // Cross product of two quaternions
   function(q) {
      return Quat(0, this.c*q.d-this.d*q.c,
                  this.d*q.b-this.b*q.d, this.b*q.c-this.c*q.b) ;
   },
   'dot': // dot product of two quaternions
   function(q) {
      return this.b*q.b+this.c*q.c+this.d*q.d ;
   },
   'normalize': // make the magnitude be 1
   function() {
      var d = Math.sqrt(this.dot(this)) ;
      return Quat(this.a/d, this.b/d, this.c/d, this.d/d) ;
   },
   'makenormal': // make a normal vector from a plane or quat
   function() {
      return Quat(0, this.b, this.c, this.d).normalize() ;
   },
   'normalizeplane': // normalize a plane
   function() {
      var d = Math.hypot(this.b, this.c, this.d) ;
      return Quat(this.a/d, this.b/d, this.c/d, this.d/d) ;
   },
   'smul': // scalar multiplication
   function(m) {
      return Quat(this.a*m, this.b*m, this.c*m, this.d*m) ;
   },
   'sum': // quaternion sum
   function(q) {
      return Quat(this.a+q.a, this.b+q.b, this.c+q.c, this.d+q.d) ;
   },
   'angle': // quaternion angle
   function() {
      return 2 * Math.acos(this.a) ;
   },
   'invrot': // quaternion inverse rotation
   function() {
      return Quat(this.a, -this.b, -this.c, -this.d) ;
   },
   'det3x3': // calculate a 3x3 determinant
   function(a00, a01, a02, a10, a11, a12, a20, a21, a22) {
      return a00 * (a11 * a22 - a12 * a21) +
             a01 * (a12 * a20 - a10 * a22) +
             a02 * (a10 * a21 - a11 * a20) ;
   },
   'rotateplane': // rotate a plane using a quaternion
   function(q) {
      var t = q.mul(Quat(0, this.b, this.c, this.d)).mul(q.invrot()) ;
      t.a = this.a ;
      return t ;
   },
   'intersect3': // find the intersection of three planes if there is one
   function(p2, p3) {
      var det = this.det3x3(this.b, this.c, this.d,
                            p2.b, p2.c, p2.d,
                            p3.b, p3.c, p3.d) ;
      if (Math.abs(det) < eps)
         return false ;
      return Quat(0,
                  this.det3x3(this.a, this.c, this.d,
                              p2.a, p2.c, p2.d, p3.a, p3.c, p3.d)/det,
                  this.det3x3(this.b, this.a, this.d,
                              p2.b, p2.a, p2.d, p3.b, p3.a, p3.d)/det,
                  this.det3x3(this.b, this.c, this.a,
                              p2.b, p2.c, p2.a, p3.b, p3.c, p3.a)/det) ;
   },
   'solvethreeplanes': // find intersection of three planes but only if interior
   // Takes three indices into a plane array, and returns the point at the
   // intersection of all three, but only if it is internal to all planes.
   function(p1, p2, p3, planes) {
      var p = planes[p1].intersect3(planes[p2], planes[p3]) ;
      if (!p)
         return p ;
      for (var i=0; i<planes.length; i++) {
         if (i != p1 && i != p2 && i != p3) {
            var dt = planes[i].b * p.b + planes[i].c * p.c + planes[i].d * p.d ;
            if ((planes[i].a > 0 && dt > planes[i].a) ||
                (planes[i].a < 0 && dt < planes[i].a))
               return false ;
         }
      }
      return p ;
   },
   'side': // is this point close to the origin, or on one or the other side?
   function(x) {
      if (x > eps)
         return 1 ;
      if (x < -eps)
         return -1 ;
      return 0 ;
   },
   'cutfaces': // Cut a set of faces by a plane and return new set
   function(faces) {
      var that = this ; // welcome to Javascript
      var d = this.a ;
      var nfaces = [] ;
      for (var j=0; j<faces.length; j++) {
         var face = faces[j] ;
         var inout = face.map(function(_){ return that.side(_.dot(that)-d)}) ;
         var seen = 0 ;
         for (var i=0; i<inout.length; i++) {
            seen |= 1<<(inout[i]+1) ;
         }
         if ((seen & 5) == 5) { // saw both sides
            for (var s=-1; s<=1; s += 2) {
               var nface = [] ;
               for (var k=0; k<face.length; k++) {
                  if (inout[k] == s || inout[k] == 0) {
                     nface.push(face[k]) ;
                  }
                  var kk = (k + 1) % face.length ;
                  if (inout[k] + inout[kk] == 0 && inout[k] != 0) {
                     var vk = face[k].dot(this) - d ;
                     var vkk = face[kk].dot(this) - d ;
                     var r = vk / (vk - vkk) ;
                     var pt = face[k].smul(1-r).sum(face[kk].smul(r)) ;
                     nface.push(pt) ;
                  }
               }
               nfaces.push(nface) ;
            }
         } else { // no split
            nfaces.push(face) ;
         }
      }
      return nfaces ;
   },
   'faceside': // which side of a plane is a face on?
   function(p, f) {
      var d = p[0] ;
      for (var i=0; i<f.length; i++) {
         var s = this.side(face[i].dot(p)-d) ;
         if (s != 0)
            return s ;
      }
      throw "Could not determine side of plane in faceside" ;
   },
   'expandfaces': // given a set of faces, expand it by a rotation set
   function(rots, faces) {
      var nfaces = [] ;
      for (var i=0; i<rots.length; i++) {
         for (var k=0; k<faces.length; k++) {
            var face = faces[k] ;
            var nface = [] ;
            for (var j=0; j<face.length; j++) {
               q = face[j].rotateplane(rots[i]) ;
               nface.push(q) ;
            }
            nfaces.push(nface) ;
         }
      }
      return nfaces ;
   },
   'sameplane': // are two planes the same?
   function(p) {
      var a = this.normalize() ;
      var b = p.normalize() ;
      return a.dist(b) < eps || a.dist(b.smul(-1)) < eps ;
   },
   'centermassface': // calculate a center of a face by averaging points
   function(face) {
      var s = Quat(0, 0, 0, 0) ;
      for (var i=0; i<face.length; i++)
         s = s.sum(face[i]) ;
      return s.smul(1.0/face.length) ;
   },
   'makecut': // make a cut from a normal vector
   function(r) {
      return Quat(r, this.b, this.c, this.d) ;
   }
} ;

// Next we define a class that yields quaternion generators for each of
// the five platonic solids.  The quaternion generators chosen are
// chosen specifically so that the first quaternion doubles as a plane
// description that yields the given Platonic solid (so for instance, the
// cubical group and octahedral group are identical in math, but we
// give distinct representations choosing the first quaternion so that
// we get the desired figure.)  Our convention is one vertex of the
// shape points precisely down.

function PlatonicGenerator() {
   if (this instanceof PlatonicGenerator) {
      return this ;
   } else {
      return new PlatonicGenerator() ;
   }
}
PlatonicGenerator.prototype = {
   'cube':
   function() {
      var s5 = Math.sqrt(0.5) ;
      return [Quat(s5, s5, 0, 0), Quat(s5, 0, s5, 0)] ;
   },
   'tetrahedron':
   function() {
      return [Quat(0.5, 0.5, 0.5, 0.5), Quat(0.5, 0.5, 0.5, -0.5)] ;
   },
   'dodecahedron':
   function() {
      var d36 = 2 * Math.PI / 10 ;
      var dx = 0.5 + 0.3 * Math.sqrt(5) ;
      var dy = 0.5 + 0.1 * Math.sqrt(5) ;
      var dd = Math.sqrt(dx*dx+dy*dy) ;
      dx /= dd ;
      dy /= dd ;
      return [Quat(Math.cos(d36), dx*Math.sin(d36), dy*Math.sin(d36), 0),
              Quat(0.5, 0.5, 0.5, 0.5)] ;
   },
   'icosahedron':
   function() {
      var dx = 1/6 + Math.sqrt(5)/6 ;
      var dy = 2/3 + Math.sqrt(5)/3 ;
      var dd = Math.sqrt(dx*dx+dy*dy) ;
      dx /= dd ;
      dy /= dd ;
      var ang = 2 * Math.PI / 6 ;
      return [Quat(Math.cos(ang), dx*Math.sin(ang), dy*Math.sin(ang), 0),
              Quat(Math.cos(ang), -dx*Math.sin(ang), dy*Math.sin(ang), 0)] ;
   },
   'octahedron':
   function() {
      var s5 = Math.sqrt(0.5) ;
      return [Quat(0.5, 0.5, 0.5, 0.5), Quat(s5, 0, 0, s5)] ;
   },
   'closure': // compute the closure of a set of generators
   // This is quadratic in the result size.  Also, it has no protection
   // against you providing a bogus set of generators that would generate
   // an infinite group.
   function(g) {
      var q = [Quat(1, 0, 0, 0)] ;
      for (var i=0; i<q.length; i++) {
         for (var j=0; j<g.length; j++) {
            var ns = q[i].mul(g[j]) ;
            var negns = ns.smul(-1) ;
            var seen = false ;
            for (var k=0; k<q.length; k++) {
               if (ns.dist(q[k]) < eps ||
                   negns.dist(q[k]) < eps) {
                  seen = true ;
                  break ;
               }
            }
            if (!seen) {
               q.push(ns) ;
            }
         }
      }
      return q ;
   },
   'uniqueplanes': // compute unique plane rotations
   // given a rotation group and a plane, find the rotations that
   // generate unique planes.  This is quadratic in the return size.
   function(p, g) {
      var planes = [] ;
      var planerot = [] ;
      for (var i=0; i<g.length; i++) {
         var p2 = p.rotateplane(g[i]) ;
         var seen = false ;
         for (var j=0; j<planes.length; j++) {
            if (p2.dist(planes[j]) < eps) {
               seen = true ;
               break ;
            }
         }
         if (!seen) {
            planes.push(p2) ;
            planerot.push(g[i]) ;
         }
      }
      return planerot ;
   },
   'getface': // compute a face given a set of planes
   // The face returned will be a set of points that lie in the first plane
   // in the given array, that are on the surface of the polytope defined
   // by all the planes, and will be returned in clockwise order.
   // This is O(planes^2 * return size + return_size^2).
   function(planes) {
      var face = [] ;
      for (var i=1; i<planes.length; i++) {
         for (var j=i+1; j<planes.length; j++) {
            var p = planes[0].solvethreeplanes(0, i, j, planes) ;
            if (p) {
               var seen = false ;
               for (var k=0; k<face.length; k++) {
                  if (p.dist(face[k]) < eps) {
                     seen = true ;
                     break ;
                  }
               }
               if (!seen) {
                  face.push(p) ;
               }
            }
         }
      }
      while (true) {
         var changed = false ;
         for (var i=0; i<face.length; i++) {
            var j = (i + 1) % face.length ;
            if (planes[0].dot(face[i].cross(face[j])) < 0) {
               var t = face[i] ;
               face[i] = face[j] ;
               face[j] = t ;
               changed = true ;
            }
         }
         if (!changed)
            break ;
      }
      return face ;
   }
} ;

// test some things
var pg = PlatonicGenerator() ;
var g = pg.cube() ;

// Other choices of base planes can be interesting but we stick with the
// basic Platonic solids for now.

var baseplane = g[0] ;
var rotations = pg.closure(g) ;
console.log("We see " + rotations.length + " rotations.") ;
var baseplanerot = pg.uniqueplanes(baseplane, rotations) ;
var baseplanes = baseplanerot.map(
                           function(_){ return baseplane.rotateplane(_) }) ;
console.log("We see " + baseplanes.length + " base planes.") ;
var baseface = pg.getface(baseplanes) ;
console.log("Basic face has " + baseface.length + " vertices.") ;

// make normals for our cuts.  We assume the first two faces are
// adjacent (which they are, by the way we construct them).

var facenormal = baseplanes[0].makenormal() ;
var edgenormal = baseplanes[0].sum(baseplanes[1]).makenormal() ;
var vertexnormal = baseface[0].makenormal() ;

var cutplanes = [] ;
cutplanes.push(facenormal.makecut(1/3)) ;
var moveplanes = [] ;
var faces = [baseface] ;

// expand cutplanes by rotations.  We only change one face here.

for (var c=0; c<cutplanes.length; c++) {
   for (var i=0; i<rotations.length; i++) {
      var q = cutplanes[c].rotateplane(rotations[i]) ;
      var seen = false ;
      for (var j=0; j<moveplanes.length; j++) {
         if (q.sameplane(moveplanes[j])) {
            seen = true ;
            break ;
         }
      }
      if (!seen) {
         moveplanes.push(q) ;
         faces = q.cutfaces(faces) ;
      }
   }
}

console.log("Faces is now " + faces.length) ;
console.log("Move planes is " + moveplanes.length) ;
