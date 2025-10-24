


CFLAGS += $$(pkg-config --cflags starpu-1.4) -g -pthread -O0 -Wall
LDLIBS += $$(pkg-config --libs starpu-1.4)

TARGETS = main


all: $(TARGETS)

run: $(TARGETS)
	./main 4 80 10

clean:
	rm -f $(TARGETS) *.o
