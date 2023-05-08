intrinsic TransposeMatrixGroup(G::GrpMat) -> GrpMat
    {Return the transpose}
    Gt := sub<GL(2,BaseRing(G)) | [Transpose(g):g in Generators(G)]>;
    return Gt;
end intrinsic;

intrinsic NonSurjectivePrimes(G::GrpMat) -> SeqEnum[RngIntElt]
    {Given G the adelic image of the Galois representation, compute the non-surjective primes}
        m:=Modulus(BaseRing(G));
        return [p:p in PrimeFactors(m)|#ChangeRing(G,GF(p)) ne #GL(2,GF(p))];
end intrinsic;

intrinsic ReducedLevel(G::GrpMat) -> RngIntElt
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
        if not p in NS and Valuation(m0,p) eq 1 then
            if m0 eq p and #G eq #GL(2,GF(p)) then
                return 1;
            elif #G/#ChangeRing(G,Integers(m0 div p)) eq #GL(2,GF(p)) then 
                m0:=m0 div p;
            end if;
        end if;
    end for;
    return m0;
end intrinsic;

intrinsic VectorOrder(v::ModTupRngElt) -> RngIntElt
    {Order of a vector with entries in Z/mZ}
    m:=#Parent(v[1]);
    g:=GCD([m, Integers()!v[2], Integers()!v[1]]);
    return m div g;
end intrinsic;

//This will be obsolete
intrinsic DegreesOfPoints(G::GrpMat) -> SetMulti
    {compute all degrees of points without doing Riemann--Roch}
    Gt := TransposeMatrixGroup(G);
    m := Modulus(BaseRing(Gt));
    H := sub<GL(2,Integers(m))|Gt,-Gt!1>;
    orb := Orbits(H);
    orbm := [ m ne 2 select #x div 2 else #x : x in orb | VectorOrder(x[1]) eq m];
    degrees := SequenceToMultiset(orbm); 
    return degrees;
end intrinsic;

//TODO: fix this
intrinsic PrimitiveDegreesOfPoints(G::GrpMat) -> SetMulti
    {compute all degrees of points without doing Riemann--Roch. 
    returns only the points ...}
    Gt := TransposeMatrixGroup(G);
    m := Modulus(BaseRing(Gt));
    H := sub<GL(2,Integers(m))|Gt,-Gt!1>;
    orb := Orbits(H);
    orbm := [ m ne 2 select <m,#x div 2, x[1]> else <m,#x, x[1]> : x in orb | VectorOrder(x[1]) eq m];
    degrees := SequenceToMultiset(orbm); 
    return degrees;
end intrinsic;

intrinsic FilterByLevelMapping(allpts::SeqEnum[Tup]) -> SeqEnum[Tup]
    {Filtering levels and degrees based on the level reduction theorem of BELOV.
    If j is a sporadic point of degree d in level m then it becomes a point of
    d/deg(f) in level n where f:X1(m)-->X1(n) is the natural projection map.}

    function CachedGenus(m, A)
        if m notin Keys(A) then
            A[m] := Genus(Gamma1(m));
        end if; 
        return A[m], A;
    end function;

    function easyRiemannRoch(listofpts,A)
        nonisolated := [];
        for pt in listofpts do
            m, deg := Explode(pt);
            genusGamma1, A := CachedGenus(m,A);
            if deg ge genusGamma1 + 1 then
                Append(~nonisolated, <m, deg>); //"easy" Riemann--Roch condition
            end if;
        end for;
        return nonisolated;
    end function;

    A := AssociativeArray();

    nonisolated := easyRiemannRoch(allpts,A);
    potisolated := SequenceToMultiset(allpts) diff SequenceToMultiset(nonisolated);



    remove := {* *};

    //TODO:rewrite this part
    // for x in potisolated do //<m, deg, P> a point of degree deg on X1(m)
    //     for y in nonisolated do //<m/b, degy, Q>
    //         //H := sub<GL(2,Integers(m))|Gt,-Gt!1>; could check that we belong to the same orbit of H to fix this
    //         if IsDivisibleBy(x[1],y[1]) then
    //             b:=x[1] div y[1]; //how much you are reducing the level by
    //             deg:=b^2*&*[Rationals() | 1-1/p^2 : p in PrimeDivisors(x[1]) | p notin PrimeDivisors(y[1])];
    //             if deg eq x[2] div y[2] then Include(~remove,x); end if;
    //         end if;
    //     end for;
    // end for;
    //See https://arxiv.org/pdf/1808.04520.pdf Prop 2.2

    potisolated := potisolated diff remove; //the remaining potentially isolated

    return potisolated;

end intrinsic;


intrinsic NotIsolated(j::FldRatElt, path::Assoc: a:=[]) -> List
    {main function to check if a j invariant is sporadic}
    if #a eq 0 then
        a := [1,0,0,-36/(j-1728),-1/(j-1728)];
        assert jInvariant(EllipticCurve(a)) eq j;
    end if;
    E:=EllipticCurve(a);
    G,n,S:=FindOpenImage(path, E);
    m0:=ReducedLevel(G);
    
    if m0 eq 1 then
        return [*Sprint(j), true, {}*];
    end if;

    G0:=ChangeRing(G,Integers(m0));

    //generate a list of <m, deg> such that E is a non CM point on X1(m) of degree deg, and...
    allpoints := PrimitiveDegreesOfPoints(ChangeRing(G0,Integers(m))); //MAJOR CHANGE!!

    potisolated := FilterByLevelMapping(allpoints);

    if #potisolated gt 0 then
        return [*Sprint(j), false, potisolated*];
    else
        return [*Sprint(j), true, potisolated*];
    end if;

end intrinsic;
