AttachSpec("../OpenImage.spec");                                                      
path := OpenImageContext("../OpenImage/data-files");                                                               
AttachSpec("isolated.spec");


if assigned seq then
    SetColumns(0);
    SetAutoColumns(false);
    n := #seq;
    if n lt 4 then 
    	for i in [1 .. 4-n] do
    		seq := "0" cat seq;
    	end for;
    end if;
    seq := "x" cat seq;
    inputs := Split(Read(seq), "\n");
    for i in [1..#inputs] do
        jinv := inputs[i];
        jinv := eval jinv;
        jinv := Rationals()!jinv;
        output := NotIsolated(jinv, path);
        print Join([Sprint(elt) : elt in output], ":");
	end for;
    exit;
end if;