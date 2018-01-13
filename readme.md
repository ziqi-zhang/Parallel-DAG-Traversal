Our project code consists of 4 files: gen.cpp, cuda.cu, serial.cpp, openmp.cpp and openmp_optimise.cpp.

Gen.cpp is used to generate data. Use following command to compile it:
g++ gen.cpp -o gen
Use following command to generate data and save it into data.txt: 
./gen <number of nodes> <number of edges> <number of source nodes> <number of sink nodes> > data.txt

Serial.cpp is used to run our serial algorithm. Use following command to compile it:
g++ serial.cpp -o serial -O3
Use following command to run serial algorithm with input file data.txt:
./serial < data.txt

Openmp.cpp is used to run our openmp algorithm. Use following command to compile it:
g++ openmp.cpp -o openmp -O3 -fopenmp
Use following command to run openmp algorithm with input file data.txt:
./openmp < data.txt

Openmp_optimised.cpp is used to run our optimised openmp algorithm. Use following command to compile it:
g++ openmp_optimised.cpp -o openmp_optimised -O3 -fopenmp
Use following command to run openmp algorithm with input file data.txt:
./openmp_optimised < data.txt

Cuda.cu is our CUDA code. Use following command to compile it:
./make-cuda.sh
Use following command to run CUDA algorithm with input file data.txt:
./cuda < data.txt

There are 4 data files here: data0.txt, data1.txt, data2.txt and data3.txt. They are experiment dataset in our report.

