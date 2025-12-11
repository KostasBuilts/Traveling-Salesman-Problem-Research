CC = gcc
CFLAGS = -O3 -fopenmp -Wall -Wextra
LDFLAGS = -fopenmp -lm
TARGET = quantum_tsp

# List all source files
SOURCES = main.c common.c tsplib.c quantum.c optimization.c
OBJECTS = $(SOURCES:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $(TARGET) $(OBJECTS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET) test5.tsp

run: $(TARGET)
	./$(TARGET) -f att48.tsp -threads 8

test: $(TARGET)
	@echo "Creating test file..."
	@echo "NAME: test5" > test5.tsp
	@echo "TYPE: TSP" >> test5.tsp
	@echo "DIMENSION: 5" >> test5.tsp
	@echo "EDGE_WEIGHT_TYPE: EUC_2D" >> test5.tsp
	@echo "NODE_COORD_SECTION" >> test5.tsp
	@echo "1 0 0" >> test5.tsp
	@echo "2 1 0" >> test5.tsp
	@echo "3 2 0" >> test5.tsp
	@echo "4 1 1" >> test5.tsp
	@echo "5 0 1" >> test5.tsp
	@echo "EOF" >> test5.tsp
	./$(TARGET) -f test5.tsp

.PHONY: all clean run test