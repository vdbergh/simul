all: simul

simul: simul.c
	gcc -Wall  -O3 simul.c -lm -lpthread -lgsl -o simul

format: 
	clang-format -i simul.c

clean:
	rm simul
