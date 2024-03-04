// Code verifying 9.0.2 where it is proved that the points on X_1(17) and X_1(24) found by the algorithm are not P^1-isolated

// Copyright (c) 2024 Isolated Points Team
// This is available under the terms of the CC-BY-4.0 License.
// A copy of the CC-BY-4.0 License should be included with this program, but if
// not then see <https://creativecommons.org/licenses/by/4.0/>.

coords_jK := function(X,j,jinv,K);
    PP := Codomain(j);
    pt := PP ! [jinv];
    base_scheme := BaseScheme(j);
    pullback_scheme := Pullback(j,pt);
    difference := Difference(pullback_scheme, base_scheme);
    differenceK := BaseChange(difference,K);
    points := Points(differenceK);
    return points;
end function;


N:=17;
j:= -882216989/131072;
//First we need to find the quartic field over which a curve with X_1(17) has quartic points corresponding to our j-invariant
E := WeierstrassModel(EllipticCurve([1, 1, 0, -660, -7600]));
assert jInvariant(E) eq j;
//E is a representative of the isomorphism class of curves with this j-invariant
fac:=Factorization(DivisionPolynomial(E,17));
polnf:=fac[1,1];
//polnf generates our quartic field
assert IsCyclic(GaloisGroup(polnf));
K<w>:=NumberField(polnf);
assert Degree(K) eq 4;
assert IsIsomorphic(K, NumberField(fac[2,1]));
assert Degree(fac[3,1]) ge 5;
// we conclude K is the unique quartic field over which X_1(17) has points corresponding to j
// Now let's find the curve over K with this j-invariant and a point of order 17 (this is not necessary, but is a sanity check).
_<u>:=PolynomialRing(K);
x:=w;
d:=w^3 - 856035*w - 341748450;
E1:=QuadraticTwist(BaseChange(E,K),d);
assert #TorsionSubgroup(E1) eq 17;


// We now restart using a nicer model of E (the one in LMFDB) to make the computations slightly nicer
j:= -882216989/131072;
//First we need to find the quartic field over which a curve with X_1(17) has quartic points corresponding to our j-invariant
E := EllipticCurve([1, 1, 0, -660, -7600]);
assert jInvariant(E) eq j;
//E is a representative of the isomorphism class of curves with this j-invariant
fac:=Factorization(DivisionPolynomial(E,17));
polnf:=fac[1,1];
//polnf generates our quartic field
assert IsCyclic(GaloisGroup(polnf));
K<w>:=NumberField(polnf);
_<x,y>:=PolynomialRing(Rationals(),2);
// We take the equation for X_1(17) from Drew Sutherland's webpage
p:=y^4 + (x^3 + x^2 - x + 2)*y^3 + (x^3 - 3*x + 1)*y^2 - (x^4 + 2*x)*y + x^3 + x^2;
A<a,b>:=AffineSpace(Rationals(),2);
X:=Curve(A,p);
r:=(x^2 + x - y)/(x^2 + x*y + x - y^2 - y);
s:=(x + 1)/(x + y + 1);
E:=EllipticCurve([s-r*s+1,r*s-r^2*s,r*s-r^2*s,0,0]);
num:=Numerator(jInvariant(E));
denom:=Denominator(jInvariant(E));
jmap := map<X -> ProjectiveSpace(Rationals(),1) | [num, denom] >;
jinv:=j;
pts:=coords_jK(X,jmap,jinv,K);
Xk:=BaseChange(X,K);
assert #pts eq 8;
// We get 8 points, they live in 2 orbits. We need to compute the RR spaces of the sum of the points in each orbit. 
P1:=Xk![pts[1,1], pts[1,2]];
auts:=Automorphisms(K);
assert auts[1](w) eq w;
// So auts[1] is the identity;
s2:=auts[2]; s3:=auts[3]; s4:=auts[4];
P2:=Xk![s2(pts[1,1]), s2(pts[1,2])];
P3:=Xk![s3(pts[1,1]), s3(pts[1,2])];
P4:=Xk![s4(pts[1,1]), s4(pts[1,2])];
assert Dimension(RiemannRochSpace(Divisor(P1)+Divisor(P2)+Divisor(P3)+Divisor(P4))) eq 2;
orb1:={P1,P2,P3,P4};
i:=1;
repeat P5:=Xk![pts[i,1], pts[i,2]]; i:=i+1; until not (P5 in orb1);
P6:=Xk![s2(P5[1]), s2(P5[2])];
P7:=Xk![s3(P5[1]), s3(P5[2])];
P8:=Xk![s4(P5[1]), s4(P5[2])];
orb2:={P5,P6,P7,P8};
assert #(orb2 join orb1) eq 8; // Sanity check.
assert Dimension(RiemannRochSpace(Divisor(P5)+Divisor(P6)+Divisor(P7)+Divisor(P8))) eq 2;

N:=24;
_<x,y>:=PolynomialRing(Rationals(),2);
// We take the equation for X_1(24) from Drew Sutherland's webpage
p:=(x^2 - 3)*y^4 + x*(x^4 + 2*x^3 + 2*x^2 + 2*x - 3)*y^2 - x^4 - x^2;
A<a,b>:=AffineSpace(Rationals(),2);
X:=Curve(A,p);
q:=y;
t:=(4*x*y + 4*y^3)/(x*y^4 - 4*x*y^2 - x - 2*y^4 - 2*y^2);
E:=EllipticCurve([0,t^2-2*q*t-2,0,-(t^2-1)*(q*t+1)^2,0]);
num:=Numerator(jInvariant(E));
denom:=Denominator(jInvariant(E));
jmap := map<X -> ProjectiveSpace(Rationals(),1) | [num, denom] >;
jinv:=16778985534208729/81000;
R<x> := PolynomialRing(Rationals()); K<a> := NumberField(x^4 - 8*x^2 + 10);
pts:=coords_jK(X,jmap,jinv,K);
Xk:=BaseChange(X,K);
P1:=Xk![pts[1,1], pts[1,2]];
P2:=Xk![pts[2,1], pts[2,2]];
P3:=Xk![pts[3,1], pts[3,2]];
P4:=Xk![pts[4,1], pts[4,2]];
Dimension(RiemannRochSpace(Divisor(P1)+Divisor(P2)+Divisor(P3)+Divisor(P4)));
