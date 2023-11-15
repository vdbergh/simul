all: simul

simul: simul.c
	gcc -Wall  -O3 simul.c -lm -lpthread -lgsl -o simul

clean:
	rm simul
