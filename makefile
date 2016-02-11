all: fft fft_parallel fft_parallel_omp fft_parallel_mpi
.PHONY: all

fft: fft.cc
	g++ fft.cc -o fft -Ofast
	
fft_parallel: fft_parallel.cc
	g++ fft_parallel.cc -o fftpar -Ofast

fft_parallel_omp: fft_parallel_omp.cc
	g++ fft_parallel_omp.cc -o fftomp -fopenmp -Ofast

fft_parallel_mpi: fft_parallel_mpi.cc
	mpic++ fft_parallel_mpi.cc -o fftmpi -Ofast
