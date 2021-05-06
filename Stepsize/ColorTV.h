/*--------------------------------------------------------------------------*/
/*--------------------------- File ColorTV.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition and implementation of the class ColorTV, which is a stepsizes 
 * rule (SR) to be used within the subgradient method (SM) [see SubGrad.h].
 * The class conforms to the interface defined by the class Stepsize [see
 * Stepsize.h]. For more details we refer to:
 *
 * F. Barahona and R. Anbil <em>The Volume Algorithm: Producing Primal
 * Solutions with a Subgradient Method</em> Math. Prog. 87(3), 385--399,
 * 2000
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
 * Copyright &copy 2001 - 2020 by Antonio Frangioni.
 */
/*--------------------------------------------------------------------------*/
/*-------------------- DEFINITIONS & IMPLEMENTATION ------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ColorTV
 #define __ColorTV  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Stepsize.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS  ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_STP 0

/* If LOG_STP > 0, the ColorTV class produces a log of its activities on the
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
/** Definition of the class ColorTV. This class implements the target stepsize
    used in the original Volume algorithm. The method is based on classifying
    the iterations based on the obtained improvement of the function Fi(). */

class ColorTV : public Stepsize
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES ---------------------------------*/
/*---------------------------------------------------------------------------*/

/** @name Public types
    The ColorTV defines one main public types: */

 enum condition { green ,  ///< good step
		  yellow , ///< not bad, not good step
		  red      ///< bad step
                  };

 /**@} ----------------------------------------------------------------------*/
 
 /*--------------------------------------------------------------------------*/
 //public enum for algoritmic parameters
 
enum int_par_colorTV{
		intg=0,	//greentestinvl
		inty,	//yellowtestinvl
		intr,	//redtestinvl
		lastintCTV
};
enum dbl_par_colorTV{
	dblBZ=0,	//BetaZero
	lastdblCTV
};

// define and initialize here the vector of int parameters names
const std::vector< std::string > int_par_string={"intr", "inty", "intr"};

// define and initialize here the vector of double parameters names
const std::vector< std::string > dbl_par_string={
	"dblBZ"};

// define and initialize here the map for int parameters names
const std::map< std::string , ThinComputeInterface::idx_type> int_par_map={
	{"intg", intg},{"inty", inty}, {"intr", intr}};

// define and initialize here the map for double parameters names
const std::map< std::string , ThinComputeInterface::idx_type > dbl_pars_map ={
	{"dblBZ", dblBZ}};

// define and initialize here the default int parameters
const std::vector<int> dflt_int_par ={
	1,
	400,
	10
};

// define and initialize here the default double parameters
const std::vector<double> dflt_dbl_par ={	
	1
};


/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** Constructor of the class. Since the constructor of ColorTV is executed
   after the one of Stepsize, the following parameters specific for ColorTV
   have to be set:

     -# VarValue BetaZero      [1] initial value of beta

     -# Index greentestinvl [1] how many consecutive green iterations
                                are allowed before increasing beta

     -# Index yellowtestinvl [400] how many consecutive yellow iterations
                             are allowed before decreasing beta

     -# Index redtestinvl    [10] how many consecutive red iterations are
                             allowed before decreasing beta */

 ColorTV( SubGradSolver *slvr=nullptr): Stepsize( slvr)
 {
    if(!slvr)
        throw(std::invalid_argument("SubGradSolver not available"));
    else{
	      greentestinvl=dflt_int_par[intg];
	      yellowtestinvl=dflt_int_par[inty];
	      redtestinvl=dflt_int_par[intr];
	      BetaZero=dflt_dbl_par[dblBZ];
    }	 
  }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 inline void SetSTPLog( std::ostream *outs = 0 , const char lvl = 0 ) override;

/*--------------------------------------------------------------------------*/

 inline void Format( void ) override;

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
/*--------------------------------------------------------------------------*/

 inline void NewStep( void ) override;

/*--------------------------------------------------------------------------*/

 bool NeedsdkM1Gk( void ) override{ return( true ); }
 
 /*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR PARAMETERS --------------------*/
/*--------------------------------------------------------------------------*/
 void set_par(const idx_type par , const int value ) override {
	switch (par){
		case(intg):
			greentestinvl=value;
		break;
		case(inty):
			yellowtestinvl=value;
		break;
		case(intr):
			redtestinvl=value;
		break;
		default:
			Stepsize::set_par(par, value);
	}
}

 void set_par(const idx_type par , const double value) override{
	switch (par){
		case(dblBZ):
			BetaZero=value;
		break;
		default:
				Stepsize::set_par(par, value);
	}
}

 int get_int_par(const idx_type par) override {
	switch (par){
		case(intg):
			return(greentestinvl);
		break;
		case(inty):
			return(yellowtestinvl);
		break;
		case(intr):
			return(redtestinvl);
		break;
		default:
				return(Stepsize::get_int_par(par));
	}
}

 double get_dbl_par(const idx_type par) override {
	switch (par){
		case(dblBZ):
			return(BetaZero);
		break;
		default:
				return(Stepsize::get_dbl_par(par));
	}
}

 int get_dflt_int_par( const idx_type par ) override{
	        	if( ( par >= intg ) && (par<=lastintCTV) )
	        	  return( dflt_int_par[par] );
	        	 else
	        	  return(Stepsize::get_dflt_int_par(par));
 }
/*--------------------------------------------------------------------------*/

 double get_dflt_dbl_par( const idx_type par ) override{
	        	if( ( par >= dblBZ ) && (par<=lastdblCTV) )
	        	  return( dflt_dbl_par[par] );
	        	 else
              return(Stepsize::get_dflt_dbl_par(par));
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

 inline const condition coloring( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 Index greentestinvl;   ///< maximum number of consecutive green iteration
 Index yellowtestinvl;  ///< maximum number of consecutive yellow iteration
 Index redtestinvl;     ///< maximum number of consecutive red iteration

 Index lastgreeniter;   ///< the last green iteration
 Index lastyellowiter;  ///< the last yellow iteration
 Index lastrediter;     ///< the last red iteration

 VarValue BetaZero;
 VarValue LwrBnd;     ///< lower bound
 VarValue lastvalue;  ///< the previous function value

 VarValue FiLmbd;     ///< current FiLambda
 Index NrIter;

/*--------------------------------------------------------------------------*/
 
  //insert in factory
 SMSpp_insert_in_factory_h;
 

 };  // end( class ColorTV )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void ColorTV::Format( void )
{
 Stepsize::Format();

 lastgreeniter = lastyellowiter = lastrediter = 0;
 LwrBnd = - Inf<VarValue>();   // lower bound, target level are unknown
 lastvalue = FiLmbd = Inf<VarValue>(); // \f$ \bar{f}_{i-1} \f$ and \f$ f_i \f$
                                    // are unknown
 Beta = BetaZero;
 NrIter = 0;
 }

/*--------------------------------------------------------------------------*/

inline void ColorTV::NewStep( void )
{
 // get the function value - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiLmbd = ReadFkVal();

 if( FiLmbd == Inf<VarValue>() )
  throw( ( "ColorTV::GetStepsize: this should not happen" ) );

 // eventually, improve the target level  - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 UpdateTargetLevel();

 // which is it the color?  - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 coloring();
 NrIter++;

 // print the level   - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 STPLOG( 1 , std::endl << "           " << " level  = " << -FiLev
	               << " ~  beta = " << Beta << std::endl << "           "
	 );

 }  // end( ColorTV::NewStep )

/*--------------------------------------------------------------------------*/

inline void ColorTV::SetSTPLog( std::ostream *outs , const char lvl ) {
 Stepsize::SetSTPLog( outs , lvl );

 #if( LOG_STP )
  if( STPLLvl > 1 ) {
   *STPLog << std::endl << "ColorTV: ~ BetaZero = " << BetaZero
	   << " ~ Color (g , y , r) = (" << greentestinvl << " , "
	   << yellowtestinvl << " , " << redtestinvl << ")" << std::endl;
   }
 #endif
 }

/*--------------------------------------------------------------------------*/

inline bool ColorTV::UpdateTargetLevel( void )
{
 bool TLchgd = false;

 // the target level is updated after a change of the lower bound  - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VarValue LwrBnd_ = GetLowerbound();
 if( LwrBnd_ > LwrBnd )
  LwrBnd = LwrBnd_;

 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( LwrBnd == - Inf<VarValue>() )
  throw((
		    "ColorTV::UpdateTargetLevel: no lower bound is given" ) );

 // target level updating  - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( FiLev < -10 ) || ( FiLev > 10 ) ) {  // nicely away from 0
  if( FiLmbd <= FiLev + std::abs( FiLev ) * 0.05 ) {
   if( ( FiLmbd < 10 ) && ( FiLmbd > - 9.5 ) )
    FiLev= - 10;
   else {
    FiLev -= 0.025 * std::abs( FiLev  );
    FiLev = std::min( FiLev , FiLmbd - 0.05 * std::abs( FiLmbd ) );
    }
   TLchgd = true;
   }
  }
 else                                       // near 0, must be careful
  if( FiLmbd < 10 ) {
   FiLev = ( FiLmbd < - 9.5 ) ? - 10 :
            ( FiLmbd > 0 ? FiLmbd * 0.95 : FiLmbd * 1.05 );
   TLchgd = true;
   }

 // failure of the target level updating - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiLev < LwrBnd ) {  // change the level
  FiLev = LwrBnd;
  TLchgd = true;
  }

 return( TLchgd );

 }  // end( ColorTV::UpdateTargetLevel )

/*--------------------------------------------------------------------------*/

inline const ColorTV::condition ColorTV::coloring( void )
{
 Index ConColor;  // number of consecutive iteration with the same color

 // green, yellow, red - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 condition lastswing;
 if( ( -(GetdGk()) <= -1e-6 ) && ( lastvalue - FiLmbd >=
	       1e-6 * std::max( std::abs( ReadBestFiVal() ) , 1.0 ) ) )  {
  lastswing = green;
  lastgreeniter = NrIter;
  }
 else
  if( ( -(GetdGk())> -1e-6 ) && ( lastvalue > FiLmbd ) ) {
   lastswing = yellow;
   lastyellowiter = NrIter;
   }
  else {
   lastswing = red;
   lastrediter = NrIter;
   }

 // change beta, if needed - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 STPLOG( 1 , std::endl << "           " );

 switch( lastswing ) {
  case( green ):  // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ConColor = NrIter - std::max( lastyellowiter , lastrediter );
   STPLOG( 1 , " # green =  " << ConColor << " ~ ScPr = "
	       << GetdGk() <<  " ~ DeltaFi = " << lastvalue - FiLmbd );
   if( ( ConColor >= greentestinvl ) && ( Beta < 2 ) ) {
    lastgreeniter = lastyellowiter = lastrediter = NrIter;
    Beta *=  2.0;
    }
   break;

  case( yellow ):  //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ConColor = NrIter - std::max( lastgreeniter , lastrediter );
   STPLOG( 2 , " # yellow =  " << ConColor << " ~ ScPr = "
	       << GetdGk() <<  " ~ DeltaFi = " << lastvalue - FiLmbd );
   if( ConColor >= yellowtestinvl ) {
    lastgreeniter = lastyellowiter = lastrediter = NrIter;
    Beta *= 1.1;
    }
   break;

  case( red ):  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ConColor = NrIter - std::max( lastgreeniter , lastyellowiter );
   STPLOG( 1 , " # red =  " << ConColor << " ~ ScPr = " << GetdGk() );
   if( ( ConColor >= redtestinvl ) && ( Beta > VarValue( 5e-4 ) ) ) {
    lastgreeniter = lastyellowiter = lastrediter = NrIter;
    Beta *= 0.67;
    }
   break;

  }  // end( switch( lastswing ) )- - - - - - - - - - - - - - - - - - - - - -

 // update the last value - - - - - - - - - - - - - - - - - - - - - - - - - -

 lastvalue = ReadFiBar();
 return( lastswing );

 }  // end( ColorTV::coloring() )

/*--------------------------------------------------------------------------*/

};  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ColorTV.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File ColorTV.h -----------------------------*/
/*--------------------------------------------------------------------------*/
