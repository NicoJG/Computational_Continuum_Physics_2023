FLAGS = \
	-Iinclude

all:
	mkdir -p data figures 
	g++ $(FLAGS) main.cpp -o main
	./main
	python plot.py

cpp11:
	mkdir -p data figures 
	g++-11 $(FLAGS) main.cpp -o main
	./main
	python plot.py