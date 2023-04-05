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

intrinsic NotIsolated(a::SeqEnum[RngIntElt], j::MonStgElt, path::Assoc) -> List
    {main function to check if a j invariant is sporadic}
    E:=EllipticCurve(a);
    G,n,S:=FindOpenImage(path, E);
    G0:=ReducedLevel(G);
    G0t := sub<GL(2,Integers(#BaseRing(G0))) | [Transpose(g):g in Generators(G0)]>;
    k:=#BaseRing(G0t);
    good:=[];
    bad:=[];
    for b in Divisors(k) do
        if b gt 12 then //X1(11) is a rank 0 elliptic curve with non noncuspidal rational pts
            ww:=MinDegreeOfPoint(ChangeRing(G0t,Integers(b)));
            //if the min degree >= g+1 then the dimension of the Riemann-Roch
            //space associated to the point of min degree is at least 2 hence it
            //is not sporadic.
            genusGamma1bplus1 := (Genus(Gamma1(b))+1);
            if  ww ge genusGamma1bplus1 then 
                Append(~good,<b, ww>); //not isolated
            else 
            //otherwise we need to check if gonality/degree reduction 
            //arguments can be used
            //then we return all of the orbit sizes / 2 as long as degree < g+1
            wws := DegreesOfPotIsolatedPoints(ChangeRing(G0t,Integers(b)), genusGamma1bplus1);
            for ww in wws do
                Append(~bad,<b, ww>); //potentially isolated
            end for;
            end if;
        end if;
    end for;
    remove:={};

//Refuting levels and degrees based on the level reduction theorem of BELOV.
//If j is a sporadic point of degree d in level m then it becomes a point of
//d/deg(f) in level n where f:X1(m)-->X1(n) is the natural projection map.

    for x in bad do 
        for y in good do
            if IsDivisibleBy(x[1],y[1]) then
                b:=x[1] div y[1];
                
                deg:=b^2*&*[Rationals() | 1-1/p^2 : p in PrimeDivisors(x[1]) | p notin PrimeDivisors(y[1])];
                if deg eq x[2] div y[2] then Include(~remove,x); end if;
            end if;
        end for;
    end for;

    bad := SequenceToSet(bad);

    //Now we check using gonality arguments if the levels, degrees in the set
    //bad but not in remove can be handled using gonality arguments.
    if bad ne remove then 
        supbad := bad diff remove; //the remaining potentially isolated
        return [*j, a, false, supbad*]; //return j-invariant we already have
    end if;
    return [*j, a, true*];    
end intrinsic;
