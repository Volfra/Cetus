#!/bin/bash
for ((i=1; i <= 20; i++))
do
    file="test$i.c"
    
    printf "############################ COMPILING C ############################\n" 
    out="test$i.o"
	#gcc "$file" -o "$out"
	gcc -fsyntax-only "$file"
	
    printf "############################ EXECUTING Cetus ############################\n"
	printf "############################ $file ############################\n"
    ./cetus -profitable-omp=1 -parallelize-loops=4 "$file"
   
    printf "############################ END Cetus ############################\n"
done