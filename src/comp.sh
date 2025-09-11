#g++ -O2 -o SESTAR SESTAR.cpp
#awk '$1 ~ /^[HSI]/ || $2 > 1000' 2079.ms1 > fast_filter.ms1
#./SESTAR fast_filter.ms1

g++ -std=c++17 -O2 -Wall -c Envelope.cpp -o Envelope.o 
g++ -std=c++17 -O2 -Wall FindEnv.cpp Envelope.o -o find_envelope
