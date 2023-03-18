AttachSpec("/usr/people/hashimot/ModularCurves/equations/equations.spec");//change if you have 
                                                       //different path to the spec file
path := OpenImageContext("/usr/people/hashimot/ModularCurves/equations/OpenImage/data-files"); //change if you have 
                                                                            //different path to 
                                                                            //data-files.


/*
PrimitiveVectors:=function(m)
        V:=RSpace(Integers(m),2);
        W:=[V![d,a]:d in Divisors(m), a in [1..m]|GCD(a,d) eq 1];
        return W;
end function;
*/

//MinDegreeOfPoint:=function(G)
        //m:=Modulus(BaseRing(G));
        //H:=sub<GL(2,Integers(m))|G,-G!1>;
        //W:=PrimitiveVectors(m);
        //c:=Min({Integers()!(#H/(2*#Stabilizer(H,v))):v in W});
        //return c;
//end function;


/*
NotSporadic:=procedure(a);
	E:=EllipticCurve(a);
	G,n,S:=FindOpenImage(path, E);
	G0:=ReducedLevel(G);
	G0t := sub<GL(2,Integers(#BaseRing(G0))) | [Transpose(g):g in Generators(G0)]>;
	k:=#BaseRing(G0t);
	for b in Divisors(k) do
	if b gt 12 then
	   bool:=MinDegreeOfPoint(ChangeRing(G0t,Integers(b))) ge (Genus(Gamma1(b))+1); 
       print bool,b,MinDegreeOfPoint(ChangeRing(G0t,Integers(b))); 
    end if;
end for;
end procedure;
*/

/*
function form_p1_data(M);
    grid := Matrix(M,M,[<1,1,0>]);
    list:=[];
    counter:=0;
    for r in [1..M-1] do
        for c in [0..M-1] do
            if (Gcd(Gcd(r,c),M) eq 1) then
                if grid[r+1,c+1] eq 0 then
                    counter:=counter + 1;
                    grid[r+1,c+1]:=counter;
                    list := list cat [[r,c]];
                    for a in [1..M] do
                        if Gcd(a,M) eq 1 then
                            rr := (r * a) mod M;
                            cc := (c * a) mod M;
                            grid[rr+1,cc+1]:=grid[r+1,c+1];
                        end if;
                    end for;
                end if;
            else
                grid[r+1,c+1]:=-1;
            end if;
        end for;
    end for;
    list := list cat [[0,1]];
    grid[1,1]:=-1;
    for c in [1..M-1] do
        grid[1,c+1]:=#list;
    end for;

    return [RSpace(Integers(M),2)!v : v in list];
end function;
*/

/*#################################################################################################
The following program checks given an elliptic curve with rational j-invariant if it represents an 
isolated point of any degree and any level. The flow of execution is as follows:

Input curve (Cremona label or with coefficient array) -> Zywina's algorithm -> Level reduction -> compute min degree for all divisors of reduced level -> refute using genus/Reiman-Roch
arguments -> refute by pushing down to lower levels using BELOV -> Output 
(True if everything is refuted / false with left over levels/degrees, otherwise).

This program was written during COUNT at CIRM.
#################################################################################################*/





NonSurjectivePrimes:=function(G)
        m:=Modulus(BaseRing(G));
        return [p:p in PrimeFactors(m)|#ChangeRing(G,GF(p)) ne #GL(2,GF(p))];
end function;


//Level reduction from the output of Zywina's algorithm
ReducedLevel:=function(G)
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
                if not p in NS and Valuation(m0,p) eq 1 and #G/#ChangeRing(G,Integers(m0 div p)) eq #GL(2,GF(p)) then
                        m0:=m0 div p;
                        G:=ChangeRing(G,Integers(m0));
                end if;
        end for;
        return G,m0;
end function;


VectorOrder:=function(v)
    m:=#Parent(v[1]);
    g:=GCD([m, Integers()!v[2], Integers()!v[1]]);
    return m div g;
end function;


//compute min degree of a point i.e. the min degree of the field of definition
//of a point expressed as a 1 x 2  vector.
MinDegreeOfPoint:=function(G)
    m:=Modulus(BaseRing(G));
    H:=sub<GL(2,Integers(m))|G,-G!1>;
    orb:=Orbits(H);
    sorb:=Sort(orb,func<o1,o2|#o1-#o2>);
    s2:=[x: x in sorb| VectorOrder(x[1]) eq m];
    return (#s2[1]) div 2;
end function;



//main function to check if a j invariant is sporadic.
NotIsolated:=function(a, j);
	E:=EllipticCurve(a);
	G,n,S:=FindOpenImage(path, E);
	G0:=ReducedLevel(G);
	G0t := sub<GL(2,Integers(#BaseRing(G0))) | [Transpose(g):g in Generators(G0)]>;
	k:=#BaseRing(G0t);
	//k;
	good:=[];
	bad:=[];
	for b in Divisors(k) do
	if b gt 12 then
			ww:=MinDegreeOfPoint(ChangeRing(G0t,Integers(b)));
            //if the min degree >= g+1 then the dimension of the Reimann-Roch
            //space associated to the point of min degree is at least 2 hence it
            //is not sporadic.
			if  ww ge (Genus(Gamma1(b))+1) then 
				Append(~good,<b, ww>);
			else 
            //otherwise we need to check if gonality/degree reduction 
            //arguments can be used
			Append(~bad,<b, ww>);
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
        supbad := bad diff remove;
        return [*j, a, false, supbad*]; //return j-invariant we already have
    end if;
    return [*j, a, true*];    
end function;


//code snippet for parallelizing the execution of NotSporadic()

if assigned seq then
    SetColumns(0);
    SetAutoColumns(false);
    seq := eval seq;
    inputs := Split(Read("isolated_refined.txt"), "\n");
    input := inputs[seq];
    i := Index(input, ")");
    jinv := input[1..i]; //construct j-invariant
    ainv := eval  input[i+2 .. #input];
    output := NotIsolated(ainv, jinv);
    print Join([Sprint(elt) : elt in output], ":");
    exit;
end if;

