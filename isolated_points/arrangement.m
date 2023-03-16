




load "PotentiallyBadCurves.m";
pos_bad_levels_single := {};
for d in pot_false_data do
    if #d[4] eq 1 then Include(~pos_bad_levels_single, Random(d[4])); end if;
end for;

pos_bad_curves_single := [**];
for d in pos_bad_levels_single do
   bad_curves_d :={}; 
    for c in pot_false_data do 
        if c[4] eq {d} then Include(~bad_curves_d, <c[1],c[2]>); end if;
    end for;
            Append(~pos_bad_curves_single, <bad_curves_d, d>);
end for;

bad_more_places := {};
for d in pot_false_data do 
    if #d[4] ne 1 then Include(~bad_more_places, <d[1], d[2], d[4]>);
    end if;
end for;

IsFromBelow := function(a, b)
    if b[1] gt  a[1] then c := a; a:= b; b:=c; end if;
    if IsDivisibleBy(a[1],b[1]) then
			c:=a[1] div b[1];
			
			deg:=c^2*&*[Rationals() | 1-1/p^2 : p in PrimeDivisors(a[1]) | p
                notin PrimeDivisors(b[1])];
			if deg eq a[2] div b[2] then return true; else return false; end if;
		end if;
    return false;
end function;

find_min := function(a)
    min := 1000000;
    min_deg := 0;
    for d in a do 
        if d[1] le min then min:= d[1]; min_deg := d[2]; end if;
    end for;
        return <min, min_deg>;
end function;


 

arrange_bad_more_places := function(a)
    bad_more_places_dummy :=[**];
   for d in a do 
      r := d[3];
      dif := [**]; 
     while #r ne 0 do 
        r_temp :=r;
        arrising_from_min := {};
        min := find_min(r);
        for x in r_temp do 
           if IsFromBelow(x, min) then Include(~arrising_from_min,x); 
              Exclude(~r, x); 
           end if;
        end for; 
        Append(~dif, arrising_from_min);
     end while;
     Append(~bad_more_places_dummy, [*d[1], d[2], dif*]);
   end for;
       return bad_more_places_dummy;
end function;
         
bad_more_places := arrange_bad_more_places(bad_more_places);

pos_bad_levels_mult:={};

for d in bad_more_places do 
    Include(~pos_bad_levels_mult, find_min(d[3][1]));
end for;


print(" ");
print("Curves with only one <level, degree> possibly bad are in the variable
pos_bad_curves_single. The variable is a list of elements of the form <Set of
bad curves at with bad level d, d>.");
print(" ");
print("levels, degrees which are the only bad levels, degrees for a curve are in
the variable pos_bad_levels_single.");
print(" ");
print("Curves with multiple possibly bad <level, degree> are in the variable
bad_more_places. The variable is a list of elements of the form [coefficients,
j-invariant, list of sets where each set has elements <level, degree> such
that they are pushing down to the same level and degree.");
print("");
print("The variable pos_bad_levels_mult contains the smallest independent
levels,degrees for curves which can be bad at many levels and degrees.");
print("");
print("Here bad means that there is a possibility that the curve represents an
isolated point at that level of that degree.");

