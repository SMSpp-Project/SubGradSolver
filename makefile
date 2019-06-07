##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of NDOBlock / SubGradSolver                                         #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by NDOBlock / SubGradSolver, i.e., core SMS++ and           #
#   SUBGRADClass. These are *not* copied into $(SUBGRADBkINC): adding those -I       #
#   directives to the compile commands will have to done by whatever "main"  #
#   makefile is using this. Analogously, any external library and the        #
#   corresponding -L< libdirs > will have to be added to the final linking   #
#   command by  whatever "main" makefile is using this.                      #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by NDOBlock / SubGradSolver.                                            #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(libSUBGRADClINC) = the -I$( SUBGRADClass library )                     #
#           $(libSUBGRADClOBJ) = the SUBGRADClass library                            #
#           $(SUBGRADBkSDR)  = the directory where the source is                 #
#                                                                            #
#   Output: $(SUBGRADBkOBJ) = the final object(s) / library                      #
#           $(SUBGRADBkH)   = the .h files to include                            #
#           $(SUBGRADBkINC) = the -I$( source directory )                        #
#                                                                            #
#                                VERSION 2.00                                #
#                               13 - 11 - 2018                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Operations Research Group                         #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################


# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

SUBGRADBkOBJ =  $(SUBGRADBkSDR)SubGradSolver.o $(SUBGRADBkSDR)Stepsize.o $(SUBGRADBkSDR)Deflection.o  $(SUBGRADBkSDR)LuckFi.o



SUBGRADBkINC = -I$(SUBGRADBkSDR)

SUBGRADBkH   =  $(SUBGRADBkSDR)SubGradSolver.h  $(SUBGRADBkSDR)Stepsize.h $(SUBGRADBkSDR)Deflection.h $(SUBGRADBkSDR)LuckFi.h
# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(SUBGRADBkOBJ) $(SUBGRADBkSDR)*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -



$(SUBGRADBkSDR)SubGradSolver.o: $(SUBGRADBkSDR)SubGradSolver.cpp $(SUBGRADBkH) \
	$(SMS++OBJ) $(libSUBGRADClOBJ) 
	$(CC) -c $*.cpp -o $@ $(SUBGRADBkINC) $(SMS++INC) $(libSUBGRADClINC) $(SW) 

$(SUBGRADBkSDR)LuckFi.o: $(SUBGRADBkSDR)LuckFi.cpp $(SUBGRADBkSDR)LuckFi.h \
	$(SMS++OBJ)
	$(CC) -c $*.cpp -o $@ $(SUBGRADBkINC) $(SMS++INC) $(SW)

$(SUBGRADBkSDR)Stepsize.o: $(SUBGRADBkSDR)Stepsize.cpp $(SUBGRADBkSDR)Stepsize.h \
	$(SMS++OBJ)
	$(CC) -c $*.cpp -o $@ $(SUBGRADBkINC) $(SMS++INC) $(SW)

$(SUBGRADBkSDR)Deflection.o: $(SUBGRADBkSDR)Deflection.cpp $(SUBGRADBkSDR)Deflection.h \
	$(SMS++OBJ)
	$(CC) -c $*.cpp -o $@ $(SUBGRADBkINC) $(SMS++INC) $(SW)


	

########################## End of makefile ###################################
