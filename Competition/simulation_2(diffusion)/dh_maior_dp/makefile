CC = gcc
CFLAGS = -lm -fopenmp -O3
SRC = competicao2D-15-01-2025.c
OUT = pde
PYTHON_SCRIPT = plot2D.py

all: compile run process

compile:
	$(CC) $(SRC) -o pde $(CFLAGS)

run:
	./pde

process:
	python3 $(PYTHON_SCRIPT)

clean:
	rm -f pde
