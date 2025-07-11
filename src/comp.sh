g++ -O2 -o SESTAR SESTAR.cpp
awk '$1 ~ /^[HSI]/ || $2 > 1000' 2079.ms1 > fast_filter.ms1
