
//code snippet for parallelizing the execution of NotIsolated()
AttachSpec("/usr/people/hashimot/OpenImage/OpenImage.spec");                                                      
path := OpenImageContext("/usr/people/hashimot/OpenImage/data-files");                                                               

if assigned seq then
    SetColumns(0);
    SetAutoColumns(false);
    seq := eval seq;
    inputs := Split(Read("adelicgenus0curves.txt"), "\n");
    input := inputs[seq];
    jinv := eval input;
    output := NotIsolated(ainv, jinv);
    print Join([Sprint(elt) : elt in output], ":");
    exit;
end if;

