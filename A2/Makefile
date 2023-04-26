FLAGS = \
	-Iinclude \
	-lfftw3 \
	-lm \
	-fsanitize=address

all:
	mkdir -p data figures 
	g++ $(FLAGS) main.cpp -o main
	rm -f data/*.csv
	rm -f figures/*.png
	rm -f figures/*.gif
	./main
	python3 plot_spectrum.py

cpp11:
	mkdir -p data figures 
	g++-11 $(FLAGS) main.cpp -o main
	rm -f data/*.csv
	rm -f figures/*.png
	rm -f figures/*.gif
	./main
	python3 plot_spectrum.py

fftw_test: fftw_test.cpp
	g++-11 $(FLAGS) fftw_test.cpp -o fftw_test
	./fftw_test