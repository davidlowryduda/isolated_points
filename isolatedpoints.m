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

//This is obsolete
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

function CoveringDegree(m,n)
    assert(m mod n eq 0);
    b:=m div n;
    c:=(n le 2 and m gt 2) select 1/2 else 1;
    return c*b^2*&*[Rationals()|(1 - 1/p^2):p in PrimeFactors(b)|n mod p ne 0];
end function;


intrinsic PrimitiveDegreesOfPoints(G::GrpMat) -> Assoc
    {Return a multiset D such that D = [<m, [a1, d1>,  .. ],[<m, [<b1, e1>,  .. ]
    s.t. for each set [<m, [<a1, d1>,  .. ] the corresponding point on X1(m)
    maps down under the natural projection map to an isolated point of degree di point on X1(a1)}
    Gt := TransposeMatrixGroup(G);
    m := Modulus(BaseRing(Gt));
    H := sub<GL(2,Integers(m))|Gt,-Gt!1>;
    orbH := Orbits(H);

    degrees:={*CartesianProduct(Integers(),Parent({<1,1>}))|  *};

    for x in orbH do
        e:=VectorOrder(x[1]);
        degs_e:={};
        for d in Divisors(e) do
            ed:=e div d;
            orb_dx:=Orbit(H,d*x[1]);
            norb_dx := #orb_dx;
            res_fld_deg := e gt 2 and ed le 2 select #x/(2*norb_dx) else #x/norb_dx;
            if CoveringDegree(e,ed) eq res_fld_deg then
                if ed le 2 then
                    Include(~degs_e,<ed,norb_dx>);
                else
                    Include(~degs_e,<ed,norb_dx div 2>);
                end if;
            end if;
        end for;
        A:={d[1]:d in degs_e};
        B:={};
        while #A gt 0 do
            a:=Min(A);
            Include(~B,a);
            A:={b:b in A|b mod a ne 0};
        end while;
        degs_e:={d:d in degs_e|d[1] in B};
        Include(~degrees, <e,degs_e>);
    end for;
    return degrees;
end intrinsic;

intrinsic FilterByRiemannRoch(primitivepts::SetMulti) -> SeqEnum[Tup]
    {Filtering levels and degrees based on Riemann--Roch.}

    function CachedGenus(m, A)
        if m notin Keys(A) then
            A[m] := Genus(Gamma1(m));
        end if; 
        return A[m], A;
    end function;

    function easyRiemannRoch(listofpts,A)
        nonisolated := [];
        for x in primitivepts do 
            n, ptset := Explode(x);
            for pt in ptset do
                m, deg := Explode(pt);
                genusGamma1, A := CachedGenus(m,A);
            if deg ge genusGamma1 + 1 then //"easy" Riemann--Roch condition
                  Append(~nonisolated, x); 
            end if;
            end for;
        end for;
        return nonisolated;
    end function;

    A := AssociativeArray();

    nonisolated := easyRiemannRoch(primitivepts,A);
    potisolated := primitivepts diff SequenceToMultiset(nonisolated);

    return potisolated;

end intrinsic;

function CondensePoints(output)
    S := {* *};
    for x in output do
        for y in x[2] do
            m := Multiplicity(output, x) - Multiplicity(S, y);
            while m gt 0 do
                Include(~S, y);
                m := m-1;
            end while;
        end for;
    end for;
    return S;
end function;

intrinsic NotIsolated(j::FldRatElt, path::Assoc) -> List
    {main function to check if a j invariant is sporadic}

    CMjinv := [ -12288000, 54000, 0, 287496, 1728, 16581375, -3375, 8000, -32768, -884736, -884736000, -147197952000, -262537412640768000];
    require j notin CMjinv : "j is a CM j-invariant. All CM j-invariants are isolated.";

    E := EllipticCurveFromjInvariant(j);
    G,n,S := FindOpenImage(path, E);
    m0 := ReducedLevel(G);
    
    if m0 eq 1 then
        return [*Sprint(j), {}*];
    end if;

    G0:=ChangeRing(G,Integers(m0));

    primitivepts := PrimitiveDegreesOfPoints(G0); //MAJOR CHANGE!!
    potisolated := FilterByRiemannRoch(primitivepts);

    if #potisolated gt 0 then
        potisolated := CondensePoints(potisolated);
        return [*Sprint(j), potisolated*];
    else
        return [*Sprint(j), potisolated*];
    end if;

end intrinsic;


intrinsic NotIsolated(j::RngIntElt, path::Assoc) -> List
    {Coerce j into the rationals if it is integral}
    return NotIsolated(Rationals()!j, path);
end intrinsic;
