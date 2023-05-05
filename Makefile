
build: make.o func.o
	mkdir build
	mv task build

func.o: func.c
	gcc -c func.c

make.o: func.o
	mpicc -o task main.c func.o

clean:
	rm -r build
	rm func.o
