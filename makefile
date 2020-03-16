
CC = gcc
CFLAGS = -ggdb -Wall -pedantic
OBJS = util.o rearrangements.o permutations.o breakpoints.o algorithms.o approx/pr.o approx/psr.o approx/prt.o approx/pt.o approx/pst.o approx/psrt.o approx/wpr.o approx/wpt.o approx/wsr.o approx/wpsr.o approx/wst.o approx/wpst.o approx/wprt.o approx/wsrt.o approx/wpsrt.o approx/improved.o
BIN = prog

file: $(OBJS) main-file.o
	$(CC) -o $(BIN) $(CFLAGS) $(OBJS) main-file.o -lm

oneperm: $(OBJS) main-one-prm.o
	$(CC) -o $(BIN) $(CFLAGS) $(OBJS) main-one-prm.o -lm

clean:
	rm *.o approx/*.o

