##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of SubGradSolver                                                #
#                                                                            #
#   Note that $(SMS++INC) is assumed to include any -I directive             #
#   corresponding to external libraries needed by SMS++, at least to the     #
#   extent in which they are needed by the parts of SMS++ used by            #
#   SubGradSolver.                                                           #
#                                                                            #
#   Input:  $(CC)        = compiler command                                  #
#           $(SW)        = compiler options                                  #
#           $(SMS++INC)  = the -I$( core SMS++ directory )                   #
#           $(SMS++OBJ)  = the libSMS++ library itself                       #
#              $(SGSkSDR)  = the directory where the source is               #
#                                                                            #
#   Output: $(SGSkOBJ)  = the final object(s) / library                      #
#                   $(SGSkH)    = the .h files to include                    #
#                 $(SGSkINC)  = the -I$( source directory )                  #
#                                                                            #
#                             Antonio Frangioni                              #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

SGSkOBJ = $(SGSkSDR)/obj/SubGradSolver.o $(SGSkSDR)/obj/Stepsize.o \
	$(SGSkSDR)/obj/Deflection.o

SGSkINC = -I$(SGSkSDR)/include

SGSkH   = $(SGSkSDR)/include/SubGradSolver.h \
	$(SGSkSDR)/include/Stepsize.h $(SGSkSDR)/include/Deflection.h \
	$(SGSkSDR)/include/Polyak.h $(SGSkSDR)/include/ColorTV.h \
	$(SGSkSDR)/include/FumeroTV.h $(SGSkSDR)/include/Volume.h \
	$(SGSkSDR)/include/PrimalDual.h

# the C++20 module smspp.subgradsolver - - - - - - - - - - - - - - - - - - -
# [see makefile-mod of SMS++]

$(eval $(call SMSppMODrule,SGSk,$(SGSkSDR),SubGradSolver,smspp.subgradsolver,\
	$(SMS++MOD)))

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(SGSkOBJ) $(SGSkSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(SGSkSDR)/obj/SubGradSolver.o: $(SGSkSDR)/src/SubGradSolver.cpp \
	$(SGSkH) $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SGSkSDR)/src/SubGradSolver.cpp -o $@ \
	$(SGSkINC) $(SMS++INC) $(SW)

$(SGSkSDR)/obj/Stepsize.o: $(SGSkSDR)/src/Stepsize.cpp \
	$(SGSkH) $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SGSkSDR)/src/Stepsize.cpp -o $@ \
	$(SGSkINC) $(SMS++INC) $(SW)

$(SGSkSDR)/obj/Deflection.o: $(SGSkSDR)/src/Deflection.cpp \
	$(SGSkH) $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(SGSkSDR)/src/Deflection.cpp -o $@ \
	$(SGSkINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
