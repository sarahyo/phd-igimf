# Start of the Makefile
# Defining variables

#======================
# Fortran 90 compiler
# (Uncomment only one)
#======================
# GNU
 F90 = gfortran
#---------------------
# Intel
#F90 = ifort
#---------------------
 Base= $(shell cd .. && pwd)
# Makefile

out:
	$(shell cd $(Base)/src && $(F90) $(Base)/src/*.f90)
	$(shell mv $(Base)/src/a.out $(Base)/bin/IGIMF)

# Cleaning everything
clean:
#	 rm $(Base)/bin/*.data 
#	 rm $(Base)/bin/*.mrv
	 rm $(Base)/bin/*.txt
	 rm $(Base)/bin/fort.*
	 rm  $(Base)/bin/IGIMF
#	 rm $(objects)
# End of the makefile
