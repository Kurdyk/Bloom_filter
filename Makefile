all: src/hash.cpp
	g++ src/hash.cpp src/fastaProcessing.cpp -o hash

clean:
	rm hash
