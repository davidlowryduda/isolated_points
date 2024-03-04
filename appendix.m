// This code verifies the claim in the proof of the Theorem in the Appendix that the ramification points X_0(37)->J_0(37)^- have j-invariant 287496.

// Copyright (c) 2024 Isolated Points Team
// This is available under the terms of the CC-BY-4.0 License.
// A copy of the CC-BY-4.0 License should be included with this program, but if
// not then see <https://creativecommons.org/licenses/by/4.0/>.


C := SmallModularCurve(37);
G,b,_:=AutomorphismGroup(C);
for g in G do;
	Cquot, f := CurveQuotient(AutomorphismGroup(C,[b(g)]));
	if Genus(Cquot) ne 1 then continue; end if;
	//X_0(37)^+ has j invariant 110592/37
	//J_0(37)^- has j invariant 1404928000/50653
	if jInvariant(Cquot) eq 110592/37 then continue; end if;
	assert jInvariant(Cquot) eq 1404928000/50653;
	ram := Support(RamificationDivisor(f))[1];
	assert Degree(ram) eq 2;
	j := jInvariant(ram,37);
	print "j", j, MinimalPolynomial(j);
end for;
