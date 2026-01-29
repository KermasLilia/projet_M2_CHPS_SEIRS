all: seirs_c_euler seirs_c_rk4
seirs_c_euler: seirs_c_euler.c
	gcc seirs_c_euler.c -o seirs_c_euler -O3
seirs_c_rk4: seirs_c_rk4.c
	gcc seirs_c_rk4.c -o seirs_c_rk4 -O3
clean:
	rm -f seirs_c_euler seirs_c_rk4
