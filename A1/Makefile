FLAGS = \
	-Iinclude

all:
	mkdir -p data figures 
	rm -f data/*.csv
	rm -f figures/*.pdf
	g++ $(FLAGS) main.cpp -o main
	./main
	python3 plot.py

cpp11:
	mkdir -p data figures
	rm -f data/*.csv
	rm -f figures/*.pdf
	g++-11 $(FLAGS) main.cpp -o main
	./main
	python3 plot.py