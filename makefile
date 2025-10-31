objects = utils_mod.o di_tri_log_mod.o k_crossings_mod.o lewin_f_terms_mod.o dilog_terms_mod.o ln_terms_mod.o N_type_branch_cut_corrs_mod.o main.o lnln_pole_int.o lnln_pole_int_a_eq_c.o

options = -fdefault-real-8 -fdefault-integer-8 -g -O0 -fbounds-check -ffpe-trap=invalid -ffpe-trap=zero -ffpe-trap=overflow -fimplicit-none -Wall -Wcompare-reals #-Wno-unused-variable  
compiler = gfortran

eval_int: $(objects)
	$(compiler) $(options) -o eval_int $(objects)

$(objects): %.o: %.f90
	$(compiler) $(options) -c $^ -o $@

clean:
	rm *.o
	rm *.mod
	rm eval_int
