MPICC=mpicc
CFLAGS=-W
LIBS=-lm -omp

TARGETS=ivp
OBJS=ivp.o mpi_ivp.o

all: $(TARGETS)

%.o: %.c
	$(MPICC) -c $< -o $@ $(CFLAGS) $(LIBS)

$(TARGETS): $(OBJS)
	$(MPICC) $^ -o $@ $(CFLAGS) $(LIBS)

clean:
	rm -f $(TARGETS) $(OBJS)
