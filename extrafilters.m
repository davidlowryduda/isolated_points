//Special cases / more filters

// Copyright (c) 2024 Isolated Points Team
// This is available under the terms of the CC-BY-4.0 License.
// A copy of the CC-BY-4.0 License should be included with this program, but if
// not then see <https://creativecommons.org/licenses/by/4.0/>.

function mptoP1(Z);
    j:=0;
    k:=0;
    N:=[];
    for i:=1 to #Z do
        add:=true;
        B:=SetToSequence(Z[i,4]);
        m,t:=Min([a[1]: a in B]);
        E:=EllipticCurve(Z[i,1]);
        G:=FindOpenImage(E);
        G0:=ChangeRing(G,Integers(m));
        if GL2Genus(G0) eq 0 then
            Gx:=GL2Borel1(BaseRing(G0));
            G1:=sub<GL(2,BaseRing(G0))|Gx,-Gx!1>;
            G0:=sub<GL(2,BaseRing(G0))|G0,-G0!1>;
             if IsConjugateSubgroup(GL(2,BaseRing(G0)), G0, G1) then
                print "Degree of the map to P1 is:", #G0 div #G1;
                print "While the input is:", B[t];
                print "";
                j:=j+1;
                add:=false;
            else
                check_el:=false;
                S:=Conjugates(GL(2,BaseRing(G0)),G1);
                for a in S do
                    Gb:=sub<GL(2,BaseRing(G0))|Generators(a) join Generators(G0)>;
                    tr, g:=IsConjugateSubgroup(GL(2,BaseRing(G0)),Gb,G1);
                    assert tr;
                    if (#Gb div #G1) eq  B[t,2] then
                        print "We have obtained an element we can eliminate!";
                        check_el:=true;
                        print "Degree of the map to P1 is:", #Gb div #G1;
                        print "While the input is:", B[t];
                        print "";
                    add:=false;
                    end if;
                end for;
                if check_el then k:=k+1; end if;
            end if;
        end if;
    if add then N:=Append(N,Z[i]); end if;
    end for;
    j;
    k;
    return N;
end function;


function gen1cases(Z);
    N:=[];
    for i:=1 to #Z do
        add:=true;
        B:=SetToSequence(Z[i,4]);
        m,t:=Min([a[1]: a in B]);
        E:=EllipticCurve(Z[i,1]);
        G:=FindOpenImage(E);
        G0:=ChangeRing(G,Integers(m));
        if GL2Genus(G0) eq 1 then
            print "We are now dealing with the case", Z[i];
            Gx:=GL2Borel1(BaseRing(G0));
            G1:=sub<GL(2,BaseRing(G0))|Gx,-Gx!1>;
            G0:=sub<GL(2,BaseRing(G0))|G0,-G0!1>;
            S:=Conjugates(GL(2,BaseRing(G0)),G1);
                for a in S do
                    Gb:=sub<GL(2,BaseRing(G0))|Generators(a) join Generators(G0)>;
                    tr, g:=IsConjugateSubgroup(GL(2,BaseRing(G0)),Gb,G1);
                    assert tr;
                    if (#Gb div #G1) eq  B[t,2] then
                        if GL2Genus(Gb) eq 0 then
                            print "We have obtained an element we can eliminate!";
                            print "Degree of the map to P1 is:", #Gb div #G1;
                            print "While the input is:", B[t];
                            print "";
                            add:=false;
                        end if;
                    end if;
                end for;
        end if;
    if add then N:=Append(N,Z[i]); end if;
    end for;
    return N;
end function;
