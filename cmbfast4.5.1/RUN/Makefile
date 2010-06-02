# Number of HDDs-1
NUMHDD=6
# NUMHDD=0
NUMPROC=14
LOCALDATADIR=/var/tmp/uros/1
# LOCALDATADIR=./d


######################################################################

default: Data0/jl.dat.0 Data0/ujl.dat.0 Data0/eos.dat.0 \
	Data0/jlens.dat.0 Data0/cmb_par

Data0/jl.dat.0: jl.dat
	lamexec n0-$(NUMHDD) duplicate jl.dat $(NUMPROC) Data0/

Data0/ujl.dat.0: ujl.dat
	lamexec n0-$(NUMHDD) duplicate ujl.dat $(NUMPROC) Data0/

Data0/eos.dat.0: eos.dat
	lamexec n0-$(NUMHDD) duplicate eos.dat $(NUMPROC) Data0/

Data0/jlens.dat.0: jlens.dat
	lamexec n0-$(NUMHDD) duplicate jlens.dat $(NUMPROC) Data0/

jl.dat:
ujl.dat:
eos.dat:
jlens.dat:

Data0/cmb_par: cmb
	lamexec n0-$(NUMHDD)  cp -r cmb Data0/cmb_par


clean:
	lamexec n0-$(NUMHDD)  rm -vf Data0/* 

.PHONY : install
install: 
	ln -sf ../cmb
	ln -sf ../jlens
	ln -sf Data0/cmb_par
	ln -sf ../jlgen
	ln -sf ../ujlgen
	ln -sf $(LOCALDATADIR)/Data0
	lamboot -v lamhosts
	lamexec n0-$(NUMHDD) createdir $(LOCALDATADIR)/Data0

