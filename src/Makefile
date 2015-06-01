CC=gcc
CUDA=nvcc
LD_FLAGS=-lm
CC_FLAGS=-O0 -g -Wall -pedantic -Wextra
CUDA_FLAGS=-O0 -g -arch=sm_35
CUDA_LD_FLAGS=-lcudart -lcudadevrt
OBJS=main.o specden.o


default: specden.x specden_cu.x

cuda: specden_cu.x

specden.x: main.o specden.o
	$(CC) $^ $(LD_FLAGS) -o $@

specden_cu.x: main.o specden_cu.o
	$(CUDA) $(CUDA_FLAGS) $^ $(LD_FLAGS) $(CUDA_LD_FLAGS) -o $@

%.o: %.c
	$(CC) $(CC_FLAGS) $< -c

specden_cu.o: specden_cu.cu
	$(CUDA) $(CUDA_FLAGS) -c $<
	#$(CUDA) $(CUDA_FLAGS) specden_rdc.o -dlink -o $@

clean:
	rm -rf *.o *.x
