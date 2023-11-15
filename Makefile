all: simul

simul: simul.c
	gcc  -O3 simul.c -lm -lpthread -lgsl -o simul

clean:
	rm simul
