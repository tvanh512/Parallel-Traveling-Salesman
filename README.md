# Parallel-Traveling-Salesman
Sequential:
g++ -o TSP_sequential TSP_sequential.cpp
./TSP_sequential

Parallel 1
mpicxx TSP_parallel_stepping_tone.cpp -o TSP_parallel_stepping_tone
mpiexec -n 2 ./TSP_parallel_stepping_tone

Parallel 2
mpicxx TSP_common_population.cpp -o TSP_common_population
mpiexec -n 2 ./TSP_common_population
