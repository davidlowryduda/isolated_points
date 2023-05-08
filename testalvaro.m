SetProfile(true);
// AttachSpec("/Users/sachihashimoto/repos/ModularCurves/equations/equations.spec");                                                      
// path := OpenImageContext("/Users/sachihashimoto/repos/ModularCurves/equations/OpenImage/data-files");                                                                           
AttachSpec("/usr/people/hashimot/ModularCurves/equations/equations.spec");                                                      
path := OpenImageContext("/usr/people/hashimot/ModularCurves/equations/OpenImage/data-files");                                                               


AttachSpec("isolated.spec");


inputs := Split(Read("alvarodata.txt"), "\n");
for  input in inputs do
    jinv := Rationals()!(eval input);
    output := NotIsolated(jinv, path);
    print Join([Sprint(elt) : elt in output], ":");
end for;