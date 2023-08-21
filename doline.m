
//code snippet for parallelizing the execution of NotIsolated()
AttachSpec("../OpenImage.spec");                                                      
path := OpenImageContext("../OpenImage/data-files");                                                               
AttachSpec("isolated.spec");
if assigned seq then
    SetColumns(0);
    SetAutoColumns(false);
    seq := eval seq;
    inputs := Split(Read("adelicgenusgt0curves.txt"), "\n");
    input := inputs[seq];
    jinv := eval input;
    jinv := Rationals()!jinv;
    output := NotIsolated(jinv, path);
    print Join([Sprint(elt) : elt in output], ":");
    exit;
end if;

