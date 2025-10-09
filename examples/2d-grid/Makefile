CFLAGS += $$(pkg-config --cflags starpu-1.4) -g -O0
LDLIBS += $$(pkg-config --libs starpu-1.4)

TARGETS = main

all: $(TARGETS)

run: $(TARGETS)
	./main 8 3 2

clean:
	rm -f $(TARGETS) *.o
