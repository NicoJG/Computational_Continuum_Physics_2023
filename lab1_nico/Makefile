all:
	g++ main_faster.cpp -o main_faster
	./main_faster

animation:
	mkdir -p data figures 
	g++ main.cpp -o main
	./main
	python plot.py