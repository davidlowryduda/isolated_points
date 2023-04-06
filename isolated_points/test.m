AttachSpec("/usr/people/hashimot/ModularCurves/equations/equations.spec");                                                      //different path to the spec file
path := OpenImageContext("/usr/people/hashimot/ModularCurves/equations/OpenImage/data-files");                                                                           
AttachSpec("isolated.spec");

if assigned seq then
    SetColumns(0);
    SetAutoColumns(false);
    seq := eval seq;
    inputs := Split(Read("isolated_refined.txt"), "\n");
    input := inputs[seq];
    i := Index(input, ")");
    jinv := input[1..i]; //construct j-invariant
    ainv := eval input[i+2 .. #input];
    output := NotIsolated(ainv, jinv, path);
    print Join([Sprint(elt) : elt in output], ":");
    exit;
end if;

