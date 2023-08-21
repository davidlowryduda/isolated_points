SetProfile(true);
AttachSpec("../OpenImage/OpenImage.spec");                                                      
path := OpenImageContext("../OpenImage/data-files");                                                                           
Attach("../ell-adic-galois-images/groups/gl2.m");

AttachSpec("isolated.spec");

inputs := Split(Read("alvarodata.txt"), "\n");
for  input in inputs do
    jinv := Rationals()!(eval input);
    output := NotIsolated(jinv, path);
    print Join([Sprint(elt) : elt in output], ":");
end for;