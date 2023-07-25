//code for degree 4 point on X1(24)

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

N:=24;
_<x,y>:=PolynomialRing(Rationals(),2);
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
