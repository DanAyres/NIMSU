

g77=gfortran
flags= -fPIC



all: AdaptSG  

print-%  : ; @echo $* = $($*)

##############################################################################################################################################
#                                                       orthpol
##############################################################################################################################################

ORTHPOL_SDIR=./ORTHPOL/src
ORTHPOL_ODIR=./ORTHPOL/obj

ORTHPOL_OBJS = $(ORTHPOL_ODIR)/dgauss.o $(ORTHPOL_ODIR)/drecur.o $(ORTHPOL_ODIR)/dlancz.o $(ORTHPOL_ODIR)/sti_new.o $(ORTHPOL_ODIR)/qgp_new.o $(ORTHPOL_ODIR)/mcdis_new.o $(ORTHPOL_ODIR)/Probability_Distributions.o $(ORTHPOL_ODIR)/SubIntervals.o

Test_ORTHPOL_OBJS = $(ORTHPOL_ODIR)/drecur.o $(ORTHPOL_ODIR)/dgauss.o $(ORTHPOL_ODIR)/dlancz.o $(ORTHPOL_ODIR)/sti_new.o $(ORTHPOL_ODIR)/qgp_new.o $(ORTHPOL_ODIR)/mcdis_new.o $(ORTHPOL_ODIR)/Probability_Distributions.o $(ORTHPOL_ODIR)/SubIntervals.o $(ORTHPOL_ODIR)/ORTHPOL_INTERFACE.o

Test: $(Test_ORTHPOL_OBJS)
	$(g77) $(flags) -I$(ORTHPOL_ODIR) -J$(ORTHPOL_ODIR) $(Test_ORTHPOL_OBJS) -o test.out $(ORTHPOL_SDIR)/Test_Orthpol_Interface.f90
	
ORTHPOL: $(ORTHPOL_OBJS)
	f2py -c -m --f90flags=-I$(ORTHPOL_ODIR) orthpol  $(ORTHPOL_OBJS)  $(ORTHPOL_SDIR)/ORTHPOL_INTERFACE.f90
  
$(ORTHPOL_ODIR)/%.o: $(ORTHPOL_SDIR)/%.f
	$(g77) $(flags) -c -o $@  $<

$(ORTHPOL_ODIR)/%.o: $(ORTHPOL_SDIR)/%.f90
	$(g77) $(flags) -c -J$(ORTHPOL_ODIR) -o  $@ $<
	


##############################################################################################################################################
#                                                       Adaptive Sparse Grids
##############################################################################################################################################

SG_SDIR=./AdaptiveSG/src
SG_ODIR=./AdaptiveSG/obj

SG_OBJS = $(SG_ODIR)/sandia_rules.o

AdaptSG: Heap GG_Meth Sandia
	
Heap:
	f2py -c -m maxheap $(SG_SDIR)/MaxHeap.f90
    
GG_Meth:
	f2py -c -m ggmethods $(SG_SDIR)/Gerstner_Griebel_Methods.f90
    
Sandia: $(SG_OBJS)
	f2py -c -m --f90flags=-I$(SG_ODIR) sandiainterface $(SG_OBJS) $(SG_SDIR)/Sandia_Interface.f90
  
$(SG_ODIR)/%.o: $(SG_SDIR)/%.f90
	$(g77) $(flags) -c -J$(SG_ODIR) -o  $@ $<
	
	
	
	



clean:
	rm $(ORTHPOL_ODIR)/*.o $(ORTHPOL_ODIR)/*.mod *.so *.out
