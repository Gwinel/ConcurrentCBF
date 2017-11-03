CC = g++

# Uncomment one of the following to switch between debug and opt mode
OPT = -O3 -DNDEBUG -DHTM_LOCK
#OPT = -g -ggdb

CFLAGS += -fPIC -std=c++11 -fno-strict-aliasing -Wall -c -I. -I./include -I/usr/include/ -I./src/ -I../sihle/ $(OPT)

LDFLAGS+= -Wall -lpthread -lssl -lcrypto 

LIBOBJECTS = \
	./src/hashutil.o \

HEADERS = $(wildcard src/*.h)

TEST = test

all: $(TEST)

clean:
	rm -f $(TEST) */*.o test con-test bmarks

test: example/test.o $(LIBOBJECTS)
	$(CC) example/test.o $(LIBOBJECTS) $(LDFLAGS) -o $@

con-test: example/con-test.o $(LIBOBJECTS)
	$(CC) example/con-test.o $(LIBOBJECTS) $(LDFLAGS) -o $@

bmarks: example/benchmarks.o $(LIBOBJECTS)
	$(CC) example/benchmarks.o $(LIBOBJECTS) $(LDFLAGS) -o $@

%.o: %.cc ${HEADERS} Makefile
	$(CC) $(CFLAGS) $< -o $@

