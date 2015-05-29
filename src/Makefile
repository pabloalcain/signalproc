CC=gcc
LD_FLAGS=-lm
CC_FLAGS=-O3 -Wall -pedantic -Wextra
OBJS=main.o specden.o kiss_fft.o

default: specden.x

specden.x: $(OBJS)
	$(CC) $^ $(LD_FLAGS) -o $@

%.o: %.c
	$(CC) $(CC_FLAGS) $< -c

clean:
	rm -rf *.o *.x
