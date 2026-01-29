all: seirs_abm_c seirs_abm_cpp
seirs_abm_c: seirs_abm.c
	gcc seirs_abm.c -o seirs_abm_c -O3 -lm
seirs_abm_cpp: seirs_abm.cpp
	g++ seirs_abm.cpp -o seirs_abm_cpp -O3
clean:
	rm -f seirs_abm_c seirs_abm_cpp
