CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o
	$(CC) main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out *.txt
exec :
	./main.exe