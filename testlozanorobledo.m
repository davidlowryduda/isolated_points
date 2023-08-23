//test on dataset from Table 4
//https://alozano.clas.uconn.edu/wp-content/uploads/sites/490/2021/04/lozano-robledo_torsion_bounds_v13.pdf

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