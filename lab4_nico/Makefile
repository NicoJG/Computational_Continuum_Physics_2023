all:
	mkdir -p data figures
	g++ -g -fsanitize=address main.cpp -o main
	./main
	python plot.py