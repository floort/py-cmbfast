all: libcmb libjlgen
libcmb:
	f77 -O2 -fPIC -c cmbfast4.5.1/cmbflat.F
	f77 -O2 -fPIC -c cmbfast4.5.1/cmbopen.F
	f77 -O2 -fPIC -c cmbfast4.5.1/driver.F
	f77 -O2 -fPIC -c cmbfast4.5.1/lensing.f
	f77 -O2 -fPIC -c cmbfast4.5.1/subroutines.F
	f77 -O2 -fPIC -c cmbfast4.5.1/params.f
	f77 -O2 -fPIC -c cmbfast4.5.1/recfast.f
	f77 -O2 -fPIC -c cmbfast4.5.1/dverk.f
	ld -L/usr/local/lib -lg2c -shared -soname libcmb.so.1 -o libcmb.so.1.0 \
		cmbflat.o cmbopen.o driver.o lensing.o subroutines.o params.o \
		recfast.o dverk.o
libjlgen:
	f77 -O2 -fPIC -c cmbfast4.5.1/jlgen.F
	ld -L/usr/local/lib -lg2c -shared -soname libjlgen.so.1 -o libjlgen.so.1.0 \
		jlgen.o
    
clean:
	rm *.o
