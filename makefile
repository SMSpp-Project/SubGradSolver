##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of SubGradSolver                                                #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by SubGradSolver, i.e., core SMS++.                     #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by SubGradSolver.                                                   #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(SBGSLVSDR)   = the directory where the source is               #
#                                                                            #
#   Output: $(SBGSLVOBJ)   = the final object(s) / library                   #
#           $(SBGSLVH)     = the .h files to include                         #
#           $(SBGSLVINC)   = the -I$( source directory )                     #
#           $(SBGSLVDSH) = the .h files to include for the Deflection and    #
#                         Stepsize components of SubGradSolver		     #
#                                                                            #
#                                VERSION 0.1                                 #
#                               22 - 04 - 21                                 #
#                                                                            #
#                              Antonio Frangioni                             #
#                               Enrico Gorgone 				     #
#				Laura Fiorella		 	             #
#                          Operations Research Group                         #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

#CC=clang++
#CC=g++
#SW=-std=c++17



# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -
SBGSLVOBJ = $(SBGSLVSDR)SubGradSolver.o

#SBGSLVINC = -I$(SBGSLVSDR) -I$(SBGSLVSDR)Deflection -I$(SBGSLVSDR)Stepsize 

SBGSLVH = $(SBGSLVSDR)SubGradSolver.h $(SBGSLVSDR)Stepsize.h \
$(SBGSLVSDR)Deflection.h $(SBGSLVSDR)PrimalDualClass/PrimalDual.h \
$(SBGSLVSDR)Deflection/Volume.h $(SBGSLVSDR)Stepsize/ColorTV.h \
$(SBGSLVSDR)Stepsize/FumeroTV.h $(SBGSLVSDR)Stepsize/Polyak.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#clean::
#	rm -f $(SBGSLVOBJ) $(SBGSLVSDR)*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(SBGSLVSDR)SubGradSolver.o: $(SBGSLVSDR)SubGradSolver.cpp $(SBGSLVH) \
	$(SMS++OBJ)
	$(CC) -c $(SW) $*.cpp -o $@ -I$(SBGSLVSDR) -I$(SBGSLVSDR)Deflection -I$(SBGSLVSDR)Stepsize -I$(SBGSLVSDR)PrimalDualClass $(SMS++INC)

########################## End of makefile ###################################
