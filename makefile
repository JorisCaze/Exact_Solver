output: main.o functions.o
	g++ -o output main.o functions.o

main.o: main.cpp
	g++ -c -Wall -Wextra main.cpp

functions.o: functions.cpp functions.h
	g++ -c -Wall -Wextra functions.cpp

clean:
	rm *.o output
