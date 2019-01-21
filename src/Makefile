hf_run: main.o molecule.o my_math_functions.o basis_function.o hartree_fock_engine.o
	g++ main.o molecule.o my_math_functions.o basis_function.o hartree_fock_engine.o -o hf_run -llapack -lblas

my_math_functions.o: my_math_functions.cpp my_math_functions.h
	g++ -c my_math_functions.cpp

main.o: main.cpp my_matrix.h
	g++ -c main.cpp

molecule.o: molecule.cpp molecule.h
	g++ -c molecule.cpp 

basis_function.o: basis_function.cpp basis_function.h
	g++ -c basis_function.cpp

hartree_fock_engine.o: hartree_fock_engine.cpp hartree_fock_engine.h
	g++ -c hartree_fock_engine.cpp

clean:
	rm *.o hf_run 
