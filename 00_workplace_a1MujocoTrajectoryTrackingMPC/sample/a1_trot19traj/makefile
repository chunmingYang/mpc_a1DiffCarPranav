COMMON=-O2 -I../../include -L../../bin -mavx -pthread -Wl,-rpath,'$$ORIGIN'
LDFLAGS = -lgsl -lgslcblas -lm -lnlopt 
CFLAGS=-c -Wno-implicit-function-declaration
CC=gcc

ROOT = a1_trot19traj
MODEL = model
# a1_s.o a1_d.o sdlib.o a1_sduforce.o
all: a1_all

a1_all: a1_main.o a1_s.o a1_d.o sdlib.o a1_sduforce.o
	$(CC) $(COMMON) a1_main.o a1_s.o a1_d.o sdlib.o a1_sduforce.o -w -lmujoco200 -lglew -lGL $(LDFLAGS) ../../bin/libglfw.so.3 ../../bin/libmujoco200.so ../../bin/libglew.so -o ../../bin/$(ROOT)

a1_main.o: a1_main.c common/params.h common/my_controller.c
	$(CC) $(COMMON) $(CFLAGS) a1_main.c

a1_d.o: $(MODEL)/a1_d.c 
	$(CC) $(COMMON) $(CFLAGS) $(MODEL)/a1_d.c

a1_s.o: $(MODEL)/a1_s.c
	$(CC) $(COMMON) $(CFLAGS) $(MODEL)/a1_s.c

a1_sduforce.o: $(MODEL)/a1_sduforce.c
	$(CC) $(COMMON) $(CFLAGS) $(MODEL)/a1_sduforce.c

sdlib.o: $(MODEL)/sdlib.c
	$(CC) $(COMMON) $(CFLAGS) $(MODEL)/sdlib.c 

clean:
	rm a1_main.o #rm *.o
