/*--------------------------------------------------------------------------*/
/*----------------------- File PrimalDual.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The class PrimalDual implements two variants of the Primal-Dual subgradient
 * method (PDSM): the simple and weighted averages. These boils down to the
 * <em>simultaneous</em> choice of a deflection coefficient and a stepsize in
 * a SM. For more details we refer to:
 *
 * Y. Nesterov <em>Primal-dual subgradient methods for convex
 * optimization</em>. Math. Prog. 120, 221-259, 2009
 *
 * Because this class choses <em>both</rm> the deflection coefficient and the
 * stepsize, it conforms to <em>both</em> the interfaces defined by the class
 * Deflection [see Deflection.h] and Stepsize [see Stepsize.h].
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

#ifndef __PrimalDual
 #define __PrimalDual /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS  ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_VOL 0

/* If LOG_VOL > 0, the PrimalDual class produces a log of its activities on
   the ostream object and at the "level of verbosity" set with the method
   SetVOLLog() [see below]. */

#if( LOG_VOL )
 #define VOLLOG( l , x ) if( VOLLLvl > l ) *VOLLog << x
 #define VOLLOG2( l , c , x ) if( ( VOLLLvl > l ) && c ) *VOLLog << x
#else
 #define VOLLOG( l , x )
 #define VOLLOG2( l , c , x )
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PrimalDualClass.h"

#include <cmath>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Definition of the class PrimalDual. This class implements two variants of
    the Primal-Dual subgradient method (PDSM): simple and weighted averages,
    under a unified SM scheme. The method has a fixed stability center, and
    the direction \f$ d_i \f$ can be regarded as the average of all the
    previous weighted subgradients, i.e.:
    \f[
     d_i =  ( \sum_{k=1}^i \upsilon_k g_k) / \Delta_i \quad,\quad
     \Delta_i = \sum_{k=1}^i \upsilon_k
    \f]
    It is worth mentioning that the subgradient weights \f$\upsilon_k \f$
    can be non-vanishing as \f$k\rightarrow \infty\f$. The stepsize rule (SR)
    and the deflection rule (DR), respectively, are:
    \f[
    \alpha_i = \upsilon_i / \Delta_i \quad,\quad
    \nu_i = \Delta_i / \omega_i
    \f]
    where \f$\omega_i\f$ is chosen in a suitable way. */

class PrimalDual : public PrimalDualClass
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
 
/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC TYPES OF THE CLASS --------------------------*/
 enum prim_dual_int{
	 intaverage=0, 
	 lastprimdual,
 };
 
 // define and initialize here the vector of int parameters names
const std::vector< std::string > int_par_string={"intaverage"};

// define and initialize here the map for int parameters names
const std::map< std::string , ThinComputeInterface::idx_type> int_par_map={
	{"intaverage", intaverage}};

// define and initialize here the default int parameters
const std::vector<int> dflt_int_par ={
	0	//inttauIter
};

 using VarValue = Function::FunctionValue;
 using c_VarValue = Function::c_FunctionValue;

 using Vec_VarValue = Function::Vec_FunctionValue;
 
 using Index = Function::Index;
 using c_Index = Function::c_Index;
 using idx_type=ThinComputeInterface::idx_type;
/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

/** Constructor of the class, which derives from both Deflection and
   Stepsize. Since the constructor of PrimalDual is executed <em>after</em>
   the one of Stepsize, the following parameters specific for PrimalDual 
     -# int average [0] is 0=true if simple averages are used, otherwise
                            weighted averages are used. Only 0 or 1 values 
							are accepted.
 
   For this class, the general parameter LpsFct -defined in Stepsize.h-
   represents the gamma factor \a F used to set \f$\omega_i\f$. */

  PrimalDual( SubGradSolver * slvr=nullptr): PrimalDualClass( slvr )
  {
	  //default inizialization
      if(!slvr)
         throw(std::invalid_argument("SubGradSolver not available"));
      else
	      average=dflt_int_par[intaverage];
   }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   inline void SetVOLLog(std::ostream *outs = 0 , const char lvl = 0 ) override;

   inline void Format( void  ) override;

/*--------------------------------------------------------------------------*/
/*-----------------------  METHODS FOR DR and SR COMPUTATION  --------------*/
/*--------------------------------------------------------------------------*/
/** The variant requires thats DR must be computed before SR. Hence, the call
    to NewDEF comes before the one to NewStep. */

   inline void NewDEF( void ) override;

   inline void NewStep( void ) override;

/*--------------------------------------------------------------------------*/

   inline c_VarValue GetDFLCoeff( void ) override;

   inline VarValue GetStepsize( bool StepIsIncr = false ) override;
   
/*--------------------------------------------------------------------------*/
/*------------------------------ ACCESS TO PARAMETERS-----------------------*/
/*--------------------------------------------------------------------------*/
   inline void set_par_d(const idx_type par , const int value ) override;

   inline void set_par_d(const idx_type par , const double value ) override;

   inline int get_int_par_d(const idx_type par) override;

   inline double get_dbl_par_d(const idx_type par) override;

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR ALGORITHMIC PARAMEERS --------------------*/
/*--------------------------------------------------------------------------*/
   
inline int get_dflt_int_par_d( const idx_type par ) override;

inline double get_dflt_dbl_par_d( const idx_type par ) override;

/*--------------------------------------------------------------------------*/

  VarValue GetLev( void ) override
  {
   throw( std::invalid_argument( "GetLev: not allowed" ) );
   }

/*--------------------------------------------------------------------------*/

  VarValue GetBeta( void ) override
  {
   throw( std::invalid_argument( "GetBeta: not allowed" ) );
   }

/*--------------------------------------------------------------------------*/

  void SetMaxBeta( const VarValue alpha ) override
  {
   throw( std::invalid_argument( "SetMaxBeta: not allowed" ) );
   }


/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

   inline void UpdateGamma( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   VarValue delta;
   VarValue Tmpdelta;
   VarValue vk;

   VarValue StepSize;
   VarValue gamma;

   VarValue beta1;
   bool average;

   VarValue LipCnst;      // Lipschitz constant

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
  //insert in factory
	SMSpp_insert_in_factory_h;

 };  // end( class PrimalDual )
 
 using VarValue = Function::FunctionValue;
 using c_VarValue = Function::c_FunctionValue;

 using Vec_VarValue = Function::Vec_FunctionValue;
 
 using Index = Function::Index;
 using c_Index = Function::c_Index;

  using idx_type=ThinComputeInterface::idx_type;

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void PrimalDual::SetVOLLog( std::ostream *outs , const char lvl )
{
 Deflection::SetVOLLog( outs , lvl );

 #if( LOG_VOL )
  if( VOLLLvl > 1 )
   if( average )
    *VOLLog << endl << "PrimalDual: Simple averages ";
   else
    *VOLLog << endl << "PrimalDual: Weighted averages ";

  *VOLLog << "~ GammaFact = " <<  LpsFct << endl;
 #endif

} // end( PrimalDual::SetVOLLog() )
/*--------------------------------------------------------------------------*/

inline void PrimalDual::Format( void )
{
 Deflection::Format();
 VarValue LipTmp = Deflection::GetGlobalLipschitz();
 if( LipTmp != LipCnst )
  LipCnst = LipTmp;

 beta1 = 1;
 delta = 0;

 }  // end( PrimalDual::Format() )
/*--------------------------------------------------------------------------*/

inline void PrimalDual::NewStep( void )
{
 // update the Lipschitz constant - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 VarValue LipTmp =Deflection::GetGlobalLipschitz();
 if( LipTmp != LipCnst )
  LipCnst = LipTmp;

 UpdateGamma( );

 // SR computation - - -  - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 StepSize = Tmpdelta / ( gamma * beta1 );

 VOLLOG( 1, endl << "            Delta = " << Tmpdelta << " ~ Beta = "
	    << beta1 << " ~ gamma = " << gamma );

 delta = Tmpdelta;
 beta1 = ( beta1 + 1/beta1 );

 }  // end( PrimalDual::NewStep() )
/*--------------------------------------------------------------------------*/

inline VarValue PrimalDual::GetStepsize( bool StepIsIncr )
{
 if( StepIsIncr )
  throw( std::invalid_argument( "PrimalDual::GetStepsize: this is not allowed" ) );

 return( StepSize );

 }  // end( PrimalDual::GetStepsize() )
/*--------------------------------------------------------------------------*/

 inline void PrimalDual::set_par_d(const idx_type par , const int value ) {
	  switch(par){
		  case(intaverage):
			if(value==0 || value==1)
				average=value;
			else
				throw(std::invalid_argument("average can be only 0 or 1"));
			break;
		 default:
			throw(std::invalid_argument("Unknown parameter"));
	  }
  }
/*--------------------------------------------------------------------------*/
inline void PrimalDual::set_par_d(const idx_type par , const double value ) {
   throw(std::invalid_argument("Unknown parameter"));
	  }

 inline int PrimalDual::get_int_par_d(const idx_type par) {
	  switch(par){
		  case(intaverage):
			return(average);
			break;
		 default:
			throw(std::invalid_argument("Unknown parameter"));
	  }
  }

  inline double PrimalDual::get_dbl_par_d(const idx_type par) {
      throw(std::invalid_argument("Unknown parameter"));
  }
/*--------------------------------------------------------------------------*/

 inline int PrimalDual::get_dflt_int_par_d( const idx_type par ) {
	        	if( ( par >= intaverage ) && (par<=lastprimdual) )
	        	  return( dflt_int_par[par] );
	        	 else
	        	  throw(std::invalid_argument("Unknown parameter"));
 }

 inline double PrimalDual::get_dflt_dbl_par_d( const idx_type par ) {
    throw(std::invalid_argument("Unknown parameter"));
 }
/*--------------------------------------------------------------------------*/

inline void PrimalDual::NewDEF( void )
{
 if( average )
  vk = 1;
 else
  vk = 1 / std::sqrt( Deflection::GetGiNorm() );

 Tmpdelta = delta + vk;

 }  // end( PrimalDual::NewDEF() )
/*--------------------------------------------------------------------------*/

inline c_VarValue PrimalDual::GetDFLCoeff( void )
{
 return( vk / Tmpdelta );
 }
/*--------------------------------------------------------------------------*/

inline void PrimalDual::UpdateGamma( void )
{
 VarValue tStar;
 tStar=GettStar();

 if( average )
  gamma = LpsFct * ( std::sqrt( LipCnst / 2.0 )  / tStar );
 else
  gamma = LpsFct / ( std::sqrt( 2.0 * LipCnst ) * tStar );

 }  // end( PrimalDual::UpdateGamma() )
/*--------------------------------------------------------------------------*/

};  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* PrimalDual.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File PrimalDual.h -------------------------------*/
/*--------------------------------------------------------------------------*/
