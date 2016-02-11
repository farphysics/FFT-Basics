# FFT-Basics
efficient and simple FFT algorithm; parallel version using openmp; parallel MPI FFT proof of concept; python examples of real function FFT and fast sine- and cosine transforms


THIS SOFTWARE COMES WITH ABSOLUTELY NO WARRANTY
use it for whatever you like though


fft_parallel.cc -- the raw Cooley-Tukey factorisation without any actual parallelisation

fft_parallel_omp.cc -- parallel; shared memory (openmp); this is usually faster than fft.cc on multicore machines

fft_parallel_mpi.cc -- distributed memory proof of concept
