# compiling (Windows/VS)
open folder in vs
run `cmake .` with the vs dev console
close and open the new .sln
edit jpegtran project, settings to change:
* set general->configuration type to dll
* clear advanced->target file extension
build
