//code snippet for parallelizing the execution of NotIsolated()

// Copyright (c) 2024 Isolated Points Team
// This is available under the terms of the CC-BY-4.0 License.
// A copy of the CC-BY-4.0 License should be included with this program, but if
// not then see <https://creativecommons.org/licenses/by/4.0/>.

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

