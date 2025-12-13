ifeq ($(OS),Windows_NT)
	CC = x86_64-w64-mingw32-g++
else
	CC = gcc
endif
CFLAGS = -O3 -fopenmp -Wall -Wextra
LDFLAGS = -fopenmp -lm
TARGET = quantum_tsp.elf

# List all source files
SOURCES = main.c common.c tsplib.c quantum.c optimization.c
OBJECTS = $(SOURCES:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $(TARGET) $(OBJECTS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

run: $(TARGET)
	./$(TARGET) -f att48.tsp -threads 8

.PHONY: all clean run test
