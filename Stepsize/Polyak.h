/*--------------------------------------------------------------------------*/
/*--------------------------- File Polyak.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The class Polyak implements the easiest target value stepsize rule, with a
 * constant level and the factor beta. The class conforms to the interface
 * defined by the class Deflection [see Deflection.h]. For more details we
 * refer to:
 *
 *  B.T. Polyak. <em>Minimization of unsmooth functionals. Zh. Vychisl.</em>
 *  Mat. Fiz. 9(3):509--521, (1969)
 *
  * \date 21 - 04 - 2021
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Department of Informatics  \n
 *         University of Pisa \n \n
 *
 * \author Enrico Gorgone \n
 *         Graphs and Mathematical Optimization (GOM) Group \n
 *         Department of Informatics \n
 *         Free University of Brussels \n \n
 *
 * \author Laura Fiorella
 *         Department of Mathematics \n
 *         University of Pisa \n \n
 *
 * Copyright &copy 2001 - 2020 by Antonio Frangioni, Enrico Gorgone.
 */
/*--------------------------------------------------------------------------*/
/*-------------------- DEFINITIONS & IMPLEMENTATION ------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __Polyak
 #define __Polyak  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Stepsize.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS  ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_STP 0

/* If LOG_STP > 0, the Polyak class produces a log of its activities on the
   ostream object and at the "level of verbosity" set with the method
   SetSTPLog() [see below]. */

#if( LOG_STP )
 #define STPLOG( l , x ) if( STPLLvl > l ) *STPLog << x
 #define STPLOG2( l , c , x ) if( ( STPLLvl > l ) && c ) *STPLog << x
#else
 #define STPLOG( l , x )
 #define STPLOG2( l , c , x )
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Definition of the class Polyak. This class implements  a target value
    stepsize rule, whereby \f$\beta_i\f$ and \f$f^{lev}_i\f$ are constant
    (they do not depend on $i$).                                          */

class Polyak : public Stepsize
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
 
/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC TYPES OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
 enum poly_dbl{
	 dblBeta=0, 
	 lastdblP
 };
 
 // define and initialize here the vector of double parameters names
const std::vector< std::string > dbl_par_string={"dblBeta"};

// define and initialize here the map for double parameters names
const std::map< std::string , ThinComputeInterface::idx_type> dbl_par_map={
	{"dblBeta", dblBeta}};

// define and initialize here the default double parameters
const std::vector<double> dflt_dbl_par ={
	1e-2
};

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** Constructor of the class. Since the constructor of Polyak is executed
    after the one of Stepsize, the following parameters specific for the
    Polyak have to be setted:

    -# VarValue Beta [ 1e-2 ] beta value. */

 Polyak( SubGradSolver *slvr=nullptr) : Stepsize( slvr)
 {
	if(!slvr)
    	throw(std::invalid_argument("SubGradSolver not available"));
	else
		Beta=dflt_dbl_par[dblBeta];
  }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 inline void SetSTPLog( std::ostream *outs = 0 , const char lvl = 0 ) override;

/*--------------------------------------------------------------------------*/

 void Format( void ) override{
  MaxBeta = Inf<VarValue>();
  FiLev = - Inf<VarValue>();
  }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
/*--------------------------------------------------------------------------*/

 inline void NewStep( void ) override;
 
 /*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR PARAMETERS --------------------*/
/*--------------------------------------------------------------------------*/
 void set_par(const idx_type par , const double value) override{
	switch(par){
		case(dblBeta):
			Beta=value;
		break;
		default:
			Stepsize::set_par(par, value);
	}
}

void set_par(const idx_type par , const int value) override{
	Stepsize::set_par(par, value);
}

 double get_dbl_par(const idx_type par) override {
	switch(par){
		case(dblBeta):
			return(Beta);
		break;
		default:
			return(Stepsize::get_dbl_par(par));
	}
}

 int get_int_par(const idx_type par) override {
	 return(Stepsize::get_int_par(par));
 }
/*--------------------------------------------------------------------------*/

 double get_dflt_dbl_par( const idx_type par ) override{
	        	if( ( par >= dblBeta ) && (par<=lastdblP) )
	        	  return( dflt_dbl_par[par] );
	        	 else
					return(Stepsize::get_dflt_dbl_par(par));
 }

  int get_dflt_int_par( const idx_type par ) override{
				return(Stepsize::get_dflt_int_par(par));
 }

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 inline bool UpdateTargetLevel( void ) override;
 

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
 
  //insert in factory
 SMSpp_insert_in_factory_h;
 

 };  // end( class Polyak )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void Polyak::SetSTPLog( std::ostream *outs , const char lvl )
{
 Stepsize::SetSTPLog( outs , lvl );

 #if( LOG_STP )
  if( STPLLvl > 1 ) {
   *STPLog << endl << "Polyak: Beta = "<< Beta << " ~ Chi = " <<  LpsFct
	   << endl;
   }
 #endif
 }

/*--------------------------------------------------------------------------*/

inline void Polyak::NewStep( void )
{
 // update the target value - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 UpdateTargetLevel();

 // some exception  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiLev == - Inf<VarValue>() )
  throw(( "Polyak::GetStepsize: no lower bound is given" ) );

 } // end ( Polyak::NewStep )

/*--------------------------------------------------------------------------*/

inline bool Polyak::UpdateTargetLevel( void )
{
 VarValue LwrBnd = GetLowerbound();
 if( LwrBnd > FiLev ) {
  FiLev = LwrBnd;
  return( true );
  }

 return( false );
 }

/*--------------------------------------------------------------------------*/

};  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Polyak.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File Polyak.h ------------------------------*/
/*--------------------------------------------------------------------------*/
