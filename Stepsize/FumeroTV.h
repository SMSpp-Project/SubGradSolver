/*--------------------------------------------------------------------------*/
/*--------------------------- File FumeroTV.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition and implementation of the class FumeroTV, which is a stepsizes
 *  rule (SR) to be used within the subgradient method (SM) [see SubGrad.h].
 * The class conforms to the interface defined by the class Stepsize
 * [see Stepsize.h]. For more details we refer to:
 *
 *  F. Fumero. <em>A modified subgradient algorithm for Lagrangean
 *  relaxation</em> Comput. Oper. Res. 28(1):33--52, (2001)
 *
 * \version 1.01
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

#ifndef __FumeroTV
 #define __FumeroTV  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Stepsize.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS  ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_STP 0

/* If LOG_STP > 0, the FumeroTV class produces a log of its activities on the
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
/** Definition of the class FumeroTV. This class implements a target value
    stepsize. At the beginning \f$ f^{lev} \f$ is set to the given lower bound
    and exponentially goes towards \f$ f^{rec} \f$. The method works in two
    phases. In the first one \f$ \beta_i \f$ is decremented, while in the
    second one \f$ \beta_i \f$ is still decremented but also incremented
    after a pre-set number of consecutive improving iterations. */

class FumeroTV : public Stepsize
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
 
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES ---------------------------------*/
/*---------------------------------------------------------------------------*/
enum int_tau{
		intr1=0,
		inteta1,
		inteta2,
		lastintTV
};
enum dbl_tau{
	dblsigmaMin=0,
	dblbeta0,
	lastdblTV
};

// define and initialize here the vector of int parameters names
const std::vector< std::string > int_par_string={"intr1", "inteta1", "inteta2"};

// define and initialize here the vector of double parameters names
const std::vector< std::string > dbl_par_string={
	"dblsigmaMin", "dblbeta0"};

// define and initialize here the map for int parameters names
const std::map< std::string , ThinComputeInterface::idx_type> int_par_map={
	{"intr1", intr1},{"inteta1", inteta1}, {"inteta2", inteta2}};

// define and initialize here the map for double parameters names
const std::map< std::string , ThinComputeInterface::idx_type > dbl_pars_map ={
	{"dblsigmaMin", dblsigmaMin}, {"dblbeta0", dblbeta0}};

// define and initialize here the default int parameters
const std::vector<int> dflt_int_par ={
	2,
	10,
	10
};

// define and initialize here the default double parameters
const std::vector<double> dflt_dbl_par ={
	1e-4,	
	1
};

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** Constructor of the class. Since the constructor of FumeroTV is executed
   after the one of Stepsize, the following parameters specific for the
   FumeroTV are default setted:

     -# VarValue sigmaMin  [ 1e-4 ] tolerance of the function \f$ sigma_r \f$

     -# Index r1        [ 2 ]    sigma_{r1} = 1/2

     -# VarValue beta0     [ 1 ]    initial value of beta

     -# Index eta1      [ 10 ]   threshold on number of failures (phase I)

     -# index eta2      [ 10 ]   threshold on number of failures (phase II)
     */

 FumeroTV( SubGradSolver *slvr =nullptr): Stepsize(slvr)
 {
	 //default inizializations
	if(!slvr)
        throw(std::invalid_argument("SubGradSolver not available"));
    else{
		r1=dflt_int_par[intr1];
		eta1=dflt_int_par[inteta1];
		eta2=dflt_int_par[inteta2];
		sigmaMin=dflt_dbl_par[dblsigmaMin];
		BetaZero=dflt_dbl_par[dblbeta0];
	}
  }


 /*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR PARAMETERS --------------------*/
/*--------------------------------------------------------------------------*/
void set_par(const idx_type par , const int value ) override{
	switch (par){
		case(intr1):
			r1=value;
		break;
		case(inteta1):
			eta1=value;
		break;
		case(inteta2):
			eta2=value;
		break;
		default:
			Stepsize::set_par(par,value);
	}
}

 void set_par(const idx_type par , const double value) override{
	switch (par){
		case(dblsigmaMin):
			sigmaMin=value;
		break;
		case(dblbeta0):
			BetaZero=value;
		break;
		default:
			Stepsize::set_par(par,value);
	}
}

 int get_int_par(const idx_type par) override{
	switch (par){
		case(intr1):
			return(r1);
		break;
		case(inteta1):
			return(eta1);
		break;
		case(inteta2):
			return(eta2);
		break;
		default:
			return(Stepsize::get_int_par(par));
	}
}

double get_dbl_par (const idx_type par) override{
	switch (par){
		case(dblsigmaMin):
			return(sigmaMin);
		break;
		case(dblbeta0):
			return(BetaZero);
		break;
		default:
			return(Stepsize::get_dbl_par(par));
	}
}

int get_dflt_int_par( const idx_type par ) override{
	        	if( ( par >= intr1 ) && (par<=lastintTV) )
	        	  return( dflt_int_par[par] );
	        	 else
	        	  return(Stepsize::get_dflt_int_par(par));
 }
/*--------------------------------------------------------------------------*/

double get_dflt_dbl_par( const idx_type par ) override{
	        	if( ( par >= dblsigmaMin ) && (par<=lastdblTV) )
	        	  return( dflt_dbl_par[par] );
	        	 else
					return(Stepsize::get_dflt_dbl_par(par));
}

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 VarValue FiLambda;   ///< FiLambda

 VarValue LwrBnd;     ///< lower bound
 VarValue FiRef;      ///< the best solution found so far

 VarValue sigma;      ///< current value of \f$ \sigma(r)\ f$
 VarValue sigmaMin;   ///< minimum threshold of \f$ \sigma(r) \f$

 VarValue BetaZero;

 Index rc;         ///< current value of r
 Index r1;
 Index r2;

 Index etac;       ///< last improvement of FiLev
 Index eta1;       ///< thresholds for max number of failures
 Index eta2;

/*--------------------------------------------------------------------------*/
 
 //insert in factory
 SMSpp_insert_in_factory_h;


/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

inline void Format( void ) override
{
 Stepsize::Format( );
 Beta = BetaZero;

 FiLambda = FiRef = Inf<VarValue>();
 // FiLambda and the reference value are not known
 LwrBnd = -Inf<VarValue>();  // no lower bound or target level is known

 etac = rc = 0;
 r2 = ceil( pow( ( - log( sigmaMin ) / 0.6933 ) , ( 1.0 / 3.26 ) ) * r1 );

 sigma = 1;

 }  // end ( Format() )

/*--------------------------------------------------------------------------*/

inline void NewStep( void ) override
{
 // get the function value   - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiLambda = ReadFkVal();

 if( FiLambda == Inf<VarValue>() )
  throw(  "GetStepsize: this should not happen" );

 // try to improve the target level - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool TLchgd = UpdateTargetLevel();

 // update target level- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( TLchgd )
  FiLev = sigma * LwrBnd + ( 1 - sigma ) * FiRef;

 // fix beta - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Beta < 1e-6 )
  Beta = 1e-6;

 // print the level - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 STPLOG( 1 , endl << "           " << " level  = " << -FiLev << " ~  beta = "
	     << Beta << endl << "           " );

 }  // end( NewStep() )

/*--------------------------------------------------------------------------*/

inline void SetSTPLog( std::ostream *outs , const char lvl ) override
{
 Stepsize::SetSTPLog( outs , lvl );

 #if( LOG_STP )
  if( STPLLvl > 1 ) {
   *STPLog << endl << "FumeroTV: ~ BetaZero = " << BetaZero
	   << " ~ (Eps0 , r1) = (" << sigma << " , " << r1 << ")"
	   << " ~  Eta = {" << eta1 << " , " << eta2 << "}" << endl;
   }
 #endif
 }  // end( SetSTPLog() )

/*--------------------------------------------------------------------------*/

inline bool UpdateTargetLevel( void ) override
{
 bool TLchgd = false;

 // the target level is updated if any changes of the lower bound occurs - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VarValue LwrBnd_ = GetLowerbound();
 if( LwrBnd_ > LwrBnd ) {
  LwrBnd = LwrBnd_;
  TLchgd = true;
  }

 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( LwrBnd == - Inf<VarValue>() )
  throw((
	           "UpdateTargetLevel: no lower bound is given" ) );
 VarValue FiRef_;
 if( FiRef == Inf<VarValue>() )  // not initialized yet
  FiRef_ = Inf<VarValue>();
 else
  FiRef_ = FiRef - ( 1e-6 * std::max( 1.0 , std::abs( FiRef ) ) );

 // beta updating- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiLambda <= FiRef_ ) {
  FiRef = FiLambda;
  if( rc >= r2 )
   Beta *= 2.0;
  etac = 0;
  TLchgd = true;
  }
 else  // FiLambda > FiRef_
  if( rc >= r2 ) {
   if( ++etac >= eta2 ) {
    Beta /= 2;
    etac = 0;
    }
   }
  else
   if( ++etac >= eta1 ) {
    Beta = Beta / (2.0 * Beta + 1.0);
    sigma = ExpFun( ++rc );
    TLchgd = true;
    etac = 0;
    }

 // failure of the target level updating - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiLev < LwrBnd ) {  // change the level
  FiLev = LwrBnd;
  TLchgd = true;
  }

 return( TLchgd );

 }  // end( UpdateTargetLevel() )

/*--------------------------------------------------------------------------*/

inline VarValue ExpFun( Index ri )
{
 return( exp( -0.6933 * pow( VarValue(ri) / VarValue( r1 ), 3.26 ) ) );

 }  // end( ExpFun() )

/*--------------------------------------------------------------------------*/
};
};  // end( namespace SMSpp_di_unipi_it )
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* FumeroTV.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File FumeroTV.h ----------------------------*/
/*--------------------------------------------------------------------------*/
