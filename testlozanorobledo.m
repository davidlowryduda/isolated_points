//test on dataset from Table 4
//https://alozano.clas.uconn.edu/wp-content/uploads/sites/490/2021/04/lozano-robledo_torsion_bounds_v13.pdf

// Copyright (c) 2024 Isolated Points Team
// This is available under the terms of the CC-BY-4.0 License.
// A copy of the CC-BY-4.0 License should be included with this program, but if
// not then see <https://creativecommons.org/licenses/by/4.0/>.

SetProfile(true);
AttachSpec("../OpenImage/OpenImage.spec");
path := OpenImageContext("../OpenImage/data-files");
Attach("../ell-adic-galois-images/groups/gl2.m");

AttachSpec("isolated.spec");

inputs := Split(Read("lozanorobledodata.txt"), "\n");
for  input in inputs do
    jinv := Rationals()!(eval input);
    output := NotIsolated(jinv, path);
    print Join([Sprint(elt) : elt in output], ":");
end for;
