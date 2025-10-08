CFLAGS += $$(pkg-config --cflags starpu-1.4) -g
LDLIBS += $$(pkg-config --libs starpu-1.4)

TARGETS = main

all: $(TARGETS)

clean:
	rm -f $(TARGETS) *.o
