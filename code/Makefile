CC = g++

CCFLAGS = -ansi -Wall -O6 -DNDEBUG -I../../lib
LDFLAGS = 
LOADLIBES = -lm        # the name used by the automatic linker rule

SRCS = parser_gr.cc parser_ss.cc timer.cc
HDRS = nodearc.h delta.h

all: serial.exe serialC.exe parallel.exe parallelC.exe

serial.exe: $(SRCS) $(HDRS) serial.c
	$(CC) $(CCFLAGS) -o serial.exe $(SRCS) serial.c $(LOADLIBES)

serialC.exe: $(SRCS) $(HDRS) serial.c
	$(CC) $(CCFLAGS) -DCHECKSUM -o serialC.exe $(SRCS) serial.c $(LOADLIBES)

parallel.exe: $(SRCS) $(HDRS) parallel.c
	$(CC) $(CCFLAGS) -o parallel.exe -fopenmp $(SRCS) parallel.c $(LOADLIBES)

parallelC.exe: $(SRCS) $(HDRS) parallel.c
	$(CC) $(CCFLAGS) -DCHECKSUM -o parallelC.exe -fopenmp $(SRCS) parallel.c $(LOADLIBES)

clean:
	rm -f *~ serial.exe serialC.exe parallel.exe parallelC.exe
