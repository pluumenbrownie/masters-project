
## COMPILE:

cd src/MCM/
g++ -std=c++11 -O3 -c *.cpp

cd ../../
g++ -std=c++11 -O3 -c includes/library.hpp
g++ -std=c++11 -O3 includes/main_routines.cpp main.cpp -include includes/library.hpp src/MCM/*.o -o MCM_Greedy.out


### OR:
# g++ -std=c++11 -O3 includes/main_routines.cpp main.cpp src/MCM/*.cpp -include includes/library.hpp -o MCM_Greedy.out