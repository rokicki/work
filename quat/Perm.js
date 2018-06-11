// Test some perm stuff.

"use strict" ;

function Perm(p_) {
   if (this instanceof Perm) {
      this.p = p_ ;
      this.n = p_.length ;
   } else {
      return new Perm(p_) ;
   }
}
Perm.prototype = {
   p:[], // The permutation itself
   n:0,  // length
   toString: // stringify
   function() {
      return 'Perm[' + this.p.join(' ') + ']' ;
   },
   mul: // multiply
   function(p2) {
      var c = Array(this.n) ;
      for (var i=0; i<this.n; i++)
         c[i] = p2.p[this.p[i]] ;
      return Perm(c) ;
   },
   inv: // inverse
   function() {
      var c = Array(this.n) ;
      for (var i=0; i<this.n; i++)
         c[this.p[i]] = i ;
      return Perm(c) ;
   },
   e: // identity
   function(n) {
      var c = Array(n) ;
      for (var i=0; i<n; i++)
         c[i] = i ;
      return Perm(c) ;
   },
   random: // random
   function(n) {
      var c = Array(n) ;
      for (var i=0; i<n; i++)
         c[i] = i ;
      for (var i=0; i<n; i++) {
         var j = i + Math.floor((n-i)*Math.random()) ;
         var t = c[i] ;
         c[i] = c[j] ;
         c[j] = t ;
      }
      return Perm(c) ;
   },
   compareTo: // comparison
   function(p2) {
      for (var i=0; i<this.n; i++)
         if (this.p[i] != p2.p[i])
            return this.p[i]-p2.p[i] ;
      return 0 ;
   },
}
//   Play around with strong generating sets.
var n = 0 ;
var g = [] ;
var e = Perm.prototype.e(n) ;
var sgs = [] ;
var Tk = [] ;
function resolve(p) {
   for (var i=p.p.length-1; i>=0; i--) {
      var j = p.p[i] ;
      if (j != i) {
         if (!sgs[i][j])
            return [false, p] ;
         p = p.mul(sgs[i][j].inv()) ;
         if (p.p[i] != i)
            throw "Element did not resolve i" ;
      }
   }
   return [true, p] ;
}
function check(k, p) {
   for (var i=k+1; i<p.p.length; i++)
      if (p.p[i] != i)
         throw "Bad call" ;
}
function knutha(k, p) {
   Tk[k].push(p) ;
   for (var i=0; i<sgs[k].length; i++) {
      if (!sgs[k][i])
         continue ;
      knuthb(k, sgs[k][i].mul(p)) ;
   }
}
function knuthb(k, p) {
   var j = p.p[k] ;
   if (!sgs[k][j]) {
      sgs[k][j] = p ;
      for (var i=0; i<Tk[k].length; i++)
         knuthb(k, p.mul(Tk[k][i])) ;
      return ;
   }
   var p2 = p.mul(sgs[k][j].inv()) ;
   if (p2.p[k] != k)
      throw "Bad fix? at " + k + " " + p2 + " from " + p ;
   var r = resolve(p2) ;
   if (r[0])
      return ;
   knutha(k-1, p2) ;
}
function getsgs() {
   sgs = [] ;
   Tk = [] ;
   for (var i=0; i<n; i++) {
      sgs.push([]) ;
      Tk.push([]) ;
      sgs[i][i] = e ;
   }
   for (var i=0; i<g.length; i++) {
      knutha(n-1, g[i]) ;
      var sz = 1 ;
      var tks = 0 ;
      for (var j=0; j<n; j++) {
         var cnt = 0 ;
         for (var k=0; k<n; k++)
            if (sgs[j][k])
               cnt++ ;
         tks += Tk[j].length ;
         sz *= cnt ;
      }
      console.log("After adding " + i + " sz is " + sz + " T " + tks) ;
   }
}
n = 100 ;
g = [Perm([1,3,2,4,6,5,7,0]), Perm([0,2,3,1,4,5,7,6])] ;
g = [Perm.prototype.random(n), Perm.prototype.random(n)] ;
e = Perm.prototype.e(n) ;
// console.log(g) ;
getsgs() ;
