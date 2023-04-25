FLAGS = \
	-Iinclude

all:
	mkdir -p data figures 
	g++ $(FLAGS) main.cpp -o main
	rm -f data/*.csv
	rm -f figures/*.png
	rm -f figures/*.gif
	./main
	python plot.py

cpp11:
	mkdir -p data figures 
	g++-11 $(FLAGS) main.cpp -o main
	rm -f data/*.csv
	rm -f figures/*.png
	rm -f figures/*.gif
	./main
	python plot.py