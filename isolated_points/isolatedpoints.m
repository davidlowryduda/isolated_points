intrinsic NonSurjectivePrimes(G::GrpMat) -> SeqEnum[RngIntElt]
    {Given G the adelic image of the Galois representation, compute the non-surjective primes}
        m:=Modulus(BaseRing(G));
        return [p:p in PrimeFactors(m)|#ChangeRing(G,GF(p)) ne #GL(2,GF(p))];
end intrinsic;

intrinsic ReducedLevel(G::GrpMat) -> GrpMat, RngIntElt
    {Level reduction from the output of Zywina's algorithm}
    m:=Modulus(BaseRing(G));
    NS:=Set(NonSurjectivePrimes(G));
    sE:={2,3} join NS; 
    m0:=&*[p^Valuation(m,p):p in sE];
    G:=ChangeRing(G,Integers(m0));
    for p in PrimeFactors(m0) do 
        while Valuation(m0,p) gt 1 and #G/#ChangeRing(G,Integers(m0 div p)) eq p^4 do
            m0:=m0 div p;
            G:=ChangeRing(G,Integers(m0));
        end while;
        if not p in NS and Valuation(m0,p) eq 1 and #G/#ChangeRing(G,Integers(m0 div p)) eq #GL(2,GF(p)) then //this #G/#ChangeRing(G,Integers(m0 div p)) eq #GL(2,GF(p)) is implied by other two conditions
            m0:=m0 div p;
            G:=ChangeRing(G,Integers(m0));
        end if;
    end for;
    return G,m0;
end intrinsic;

intrinsic VectorOrder(v::ModTupRngElt) -> RngIntElt
    {Order of a vector with entries in Z/mZ}
    m:=#Parent(v[1]);
    g:=GCD([m, Integers()!v[2], Integers()!v[1]]);
    return m div g;
end intrinsic;

intrinsic MinDegreeOfPoint(G::GrpMat) -> RngIntElt
    {compute min degree of a point i.e. the min degree of the field of definition
    of a point expressed as a 1 x 2  vector}
    m:=Modulus(BaseRing(G));
    H:=sub<GL(2,Integers(m))|G,-G!1>;
    orb:=Orbits(H);
    sorb:=Sort(orb,func<o1,o2|#o1-#o2>);
    s2:=[x: x in sorb| VectorOrder(x[1]) eq m]; //equivalently Minimum([x: x in orb| VectorOrder(x[1]) eq m]);
    return (#s2[1]) div 2;
end intrinsic;

intrinsic DegreesOfPotIsolatedPoints(G::GrpMat, g::RngIntElt) -> SeqEnum[RngIntElt]
    {compute all degrees of points, and returns the ones that are greater than g}
    m := Modulus(BaseRing(G));
    H := sub<GL(2,Integers(m))|G,-G!1>;
    orb := Orbits(H);
    orbm := [x : x in orb | VectorOrder(x[1]) eq m]; 
    orbmg := [#x div 2 : x in orbm | (#x div 2) ge g];
    degrees := SetToSequence(Set(orbmg)); //remove duplicates
    return degrees;
end intrinsic;

intrinsic DegreesOfPoints(G::GrpMat) -> SeqEnum[RngIntElt]
    {compute all degrees of points without doing Riemann-Roch}
    m := Modulus(BaseRing(G));
    H := sub<GL(2,Integers(m))|G,-G!1>;
    orb := Orbits(H);
    orbm := [#x div 2 : x in orb | VectorOrder(x[1]) eq m]; 
    degrees := SetToSequence(Set(orbm)); //remove duplicates
    return degrees;
end intrinsic;

//TODO: rename this function
intrinsic RefuteLevel(potisolated::SeqEnum[Tup], allpts::SeqEnum[Tup]) -> SeqEnum[Tup]
    {Refuting levels and degrees based on the level reduction theorem of BELOV.
    If j is a sporadic point of degree d in level m then it becomes a point of
    d/deg(f) in level n where f:X1(m)-->X1(n) is the natural projection map.}
	remove := {};
    for x in potisolated do 
        for y in allpts do
            if IsDivisibleBy(x[1],y[1]) then
                b:=x[1] div y[1];
                
                deg:=b^2*&*[Rationals() | 1-1/p^2 : p in PrimeDivisors(x[1]) | p notin PrimeDivisors(y[1])];
                if deg eq x[2] div y[2] then Include(~remove,x); end if;
            end if;
        end for;
    end for;

    return remove;

end intrinsic;

intrinsic NotIsolated(a::SeqEnum[RngIntElt], j::MonStgElt, path::Assoc) -> List
    {main function to check if a j invariant is sporadic}
    E:=EllipticCurve(a);
    G,n,S:=FindOpenImage(path, E);
    G0:=ReducedLevel(G);
    G0t := sub<GL(2,Integers(#BaseRing(G0))) | [Transpose(g):g in Generators(G0)]>;
    k:=#BaseRing(G0t);

    allpoints := []; //generate a list of <l, deg> such that E is a non CM point on X1(l) of degree deg 
    potisolated := [];
    for l in Divisors(k) do
        if l gt 12 then //X1(11) is a rank 0 elliptic curve with non noncuspidal rational pts

            listofdeg := DegreesOfPoints(ChangeRing(G0t,Integers(l)));
	genusGamma1lplus1 := Genus(Gamma1(l))+1;
            for deg in listofdeg do
                Append(~allpoints,<l, deg>);
            	if deg ge genusGamma1lplus1 then
                	Append(~potisolated, <l, deg>); //"easy" Riemann--Roch condition
            	end if;
	end for;
   
        end if;
    end for;

    remove := RefuteLevel(potisolated, allpoints);
	potisolated := SequenceToSet(potisolated);
	potisolated := potisolated diff remove; //the remaining potentially isolated

    if #potisolated gt 0 then
        return [*j, a, false, potisolated*]; 
    else
        return [*j, a, true, potisolated*]; 
    end if;

end intrinsic;
