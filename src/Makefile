CC=gcc
CUDA=nvcc
LD_FLAGS=-lm
CC_FLAGS=-O0 -g -Wall -pedantic -Wextra
CUDA_FLAGS=-O0 -arch=sm_30
CUDA_LD_FLAGS=-lcudart -lcudadevrt
OBJS=main.o specden.o


default: specden.x

cuda: specden_cu.x

specden.x: main.o specden.o
	$(CC) $^ $(LD_FLAGS) -o $@

specden_cu.x: main.o specden_cu.o
	$(CC) $^ $(LD_FLAGS) -lcudart -o $@

%.o: %.c
	$(CC) $(CC_FLAGS) $< -c

specden_cu.o: specden_cu.cu
	$(CUDA) $(CUDA_FLAGS) -rdc=true -c $< -o specden_rdc.o
	$(CUDA) $(CUDA_FLAGS) $(CUDA_LD_FLAGS) specden_rdc.o -dlink -o $@

clean:
	rm -rf *.o *.x
