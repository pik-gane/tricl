# tricl
TriCl model in C++

COMPILE with:

   cd build/default

   cmake ../../
   
   cmake --build .

PROFILE with 
   
   valgrind --tool=callgrind --callgrind-out-file=/tmp/callgrind.out ./tricl; kcachegrind /tmp/callgrind.out
