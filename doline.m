
//code snippet for parallelizing the execution of NotIsolated()
AttachSpec("/usr/people/hashimot/OpenImage/OpenImage.spec");                                                      
path := OpenImageContext("/usr/people/hashimot/OpenImage/data-files");                                                               
AttachSpec("isolated.spec");
if assigned seq then
    SetColumns(0);
    SetAutoColumns(false);
    seq := eval seq;
    inputs := Split(Read("alvarodata.txt"), "\n");
    input := inputs[seq];
    jinv := eval input;
    jinv := Rationals()!jinv;
    output := NotIsolated(jinv, path);
    print Join([Sprint(elt) : elt in output], ":");
    exit;
end if;

