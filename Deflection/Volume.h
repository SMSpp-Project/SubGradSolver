/*--------------------------------------------------------------------------*/
/*--------------------------- File Volume.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The class Volume implements the deflection rule of the modified Volume
 * algorithm, i.e., that of the "poorman's" Bundle method. For more details
 * we refer to:
 *
 * L. Bahiense, N. Maculan, C. Sagastizabal. <em>The volume algorithm
 * revisited: relation with bundle methods. </em> Math. Prog. 94, 41-69, 2002
 *
 * The class conforms to the interface defined by the class Deflection [see
 * Deflection.h].
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
 * Copyright &copy 2001 - 2020 by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*-------------------- DEFINITIONS & IMPLEMENTATION ------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __Volume
 #define __Volume  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS  ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_VOL 0

/* If LOG_VOL > 0, the Volume class produces a log of its activities on the
   ostream object and at the "level of verbosity" set with the method
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

#include "Deflection.h"
#include <algorithm>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Definition of the class Volume. This class implements the Volume
    algorithm. The method is revisited and incorporated in the SubGradSolver.
	The deflection coefficient is found solving a quadratic problem
    that involves the subgradient \f$ g_i \f$, the direction \f$ d_{i-1} \f$,
    and their respective linearization errors \f$ \sigma_i \f$  and
    \f$ \epsilon_{i-1} \f$ at the stability center \f$ \bar{\lambda}_i\f$,
    namely
    \f[
     \tau_i = 
     \arg\min \left\{ \; \nu_{i-1} \left\| \tau g_i 
                         + (1 - \tau) d_{i-1} \right\|^2/2
                         +  \tau \sigma_i(\bar{\lambda}_i)
                         + (1 - \tau)\epsilon_{i-1}(\bar{\lambda}_i)
        \;:\; \tau \in [0, 1]\; \right\}
    \f]  */

class Volume : public Deflection
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
 
 /*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES ---------------------------------*/
/*---------------------------------------------------------------------------*/
enum int_tau{
		inttauIter=0,	//[50]    $\tau_p$
		intlast
};
enum dbl_tau{
	dbltauInit=0,	//[0.05]  \f$ \tau_0 \f$
	dbltauMin,		//[1e-4]  safety threshold of \f$ \tau:~\tau_{\min}\f$
	dbltauFactor,	//[0.5]   factor \f$ \tau_{f} \f$
	dblm,		//[1e-3]  descent parameter */
	dblast
};

// define and initialize here the vector of int parameters names
const std::vector< std::string > int_par_string={"inttauIter"};

// define and initialize here the vector of double parameters names
const std::vector< std::string > dbl_par_string={
	"dbltauInit", "dbltauMin", "dbltauFactor", "dblm"};

// define and initialize here the map for int parameters names
const std::map< std::string , ThinComputeInterface::idx_type> int_par_map={
	{"inttauIter", inttauIter}};

// define and initialize here the map for double parameters names
const std::map< std::string , ThinComputeInterface::idx_type > dbl_pars_map ={
	{"dbltauInit", dbltauInit}, {"dbltauMin", dbltauMin}, 
	{"dbltauFactor", dbltauFactor}, {"dblm", dblm}};

// define and initialize here the default int parameters
const std::vector<int> dflt_int_par ={
	50	//inttauIter
};

// define and initialize here the default double parameters
const std::vector<double> dflt_dbl_par ={
	0.05,	//dbltauInit
	1e-4,	//dbltauMin
	0.5,	//dbltauFactor
	1e-3	//dblm
	
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
/** Constructor of the class. Each non-specified parameter is given a default 
	value, shown in [] below.

    -#   VarValue tauInit   [0.05]  \f$ \tau_0 \f$

    -#   VarValue tauMin    [1e-4]  safety threshold of \f$ \tau:~\tau_{\min} \f$

    -#   VarValue tauFactor [0.5]   factor \f$ \tau_{f} \f$

    -#   VarValue tauIter   [50]    $\tau_p$

    -#   VarValue m         [1e-3]  descent parameter */

 Volume( SubGradSolver * slvr =nullptr) : Deflection( slvr )
 {
	 //default inizialization
    if(!slvr)
        throw(std::invalid_argument("SubGradSolver not available"));
    else{
		tauInit=dflt_dbl_par[dbltauInit];
		tauMin=dflt_dbl_par[dbltauMin];
		tauFactor=dflt_dbl_par[dbltauFactor];
		tauIter=dflt_int_par[inttauIter];
		m=dflt_dbl_par[dblm]; 
	}
  }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 inline void SetVOLLog( std::ostream *outs = 0 , const char lvl = 0 ) override;

/*--------------------------------------------------------------------------*/

 void Format( void ) override
 {
  Deflection::Format();
  alpha = lastvalue = Inf<VarValue>();  // deflection coefficient is unknown
  tau = tauInit;
  }
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
/*-------------------- METHODS FOR DEFLECTION COMPUTATION ------------------*/
/*--------------------------------------------------------------------------*/

 inline void NewDEF( void ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 inline c_VarValue GetDFLCoeff( void ) override;

 inline const bool DoSS( void ) override{
 bool serious;
 VarValue delta= GetEpsLin();
  printf("delta1=%f e m=%f\n", delta, m);
 delta *=m * std::max( 1.0, std::abs(ReadBestFiVal()));
 printf("BestFi=%f\n", ReadBestFiVal());
printf("delta2=%f\n", delta);
 // if there is an improvement of Fi(), a SS occurs; otherwise, a NS occurs
 printf("FiVal- FVal=%f-%f=%f\n", ReadFiVal(),ReadFVal(), (ReadFiVal() - ReadFVal()));
 if( (ReadFiVal() - ReadFVal()) >= delta )
  serious = true;
 else
  serious = false;
printf("Serious=%d\n", serious);
 return( serious );

 }  // end( Volume::DoSS() )

 inline c_VarValue Delta( void ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 inline void Solve( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 double tauInit;     // initial tau
 double tauMin;      // the minimum value of tau
 double tauFactor;   // the decreasing factor
 int tauIter;     // the number of iterations after a decrease
 double m;           // descent parameter

 double tau;         // the safeguarded parameter tau
 VarValue lastvalue;   // last function value Fi()
 VarValue alpha;       // deflection coefficient

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
  //insert in factory
	SMSpp_insert_in_factory_h;

 };  // end( class Volume )
 using VarValue = Function::FunctionValue;
 using c_VarValue = Function::c_FunctionValue;

 using Vec_VarValue = Function::Vec_FunctionValue;
 
 using Index = Function::Index;
 using c_Index = Function::c_Index;
using idx_type=ThinComputeInterface::idx_type;

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

 inline void Volume::set_par_d(const idx_type par , const int value ) {
	/*default value shown in [] below.
		 -#   VarValue tauIter   [50]    $\tau_p$	*/
	switch(par){
		case(inttauIter):
			tauIter=value;
		break;
		default:
			throw(std::invalid_argument(("Unknown parameter")));
	}
}//end(Volume::set_par())
/*--------------------------------------------------------------------------*/

 inline void Volume::set_par_d(const idx_type par , const double value) {
	/*default value, shown in []
   below.

    -#   VarValue tauInit   [0.05]  \f$ \tau_0 \f$

    -#   VarValue tauMin    [1e-4]  safety threshold of \f$ \tau:~ \tau_{\min}\f$

    -#   VarValue tauFactor [0.5]   factor \f$ \tau_{f} \f$

    -#   VarValue m         [1e-3]  descent parameter */
	
	switch(par){
		case(dbltauInit):
			tauInit=value;
		break;
		case(dbltauMin):
			tauMin=value;
		break;
		case(dbltauFactor):
			tauFactor=value;
		break;
		case(dblm):
			m=value;
		break;
		default:
			throw(std::invalid_argument("Unknown parameter"));	
	}
}//end(Volume::set_par())
/*--------------------------------------------------------------------------*/

 inline int Volume::get_int_par_d(const idx_type par) {
	switch(par){
		case(inttauIter):
			return(tauIter);
		break;
		default:
			throw(std::invalid_argument("Unknown parameter"));
	}
}//end(Volume::get_int_par())
/*--------------------------------------------------------------------------*/

 inline double Volume::get_dbl_par_d(const idx_type par) {
	switch(par){
		case(dbltauInit):
			return(tauInit);
		break;
		case(dbltauMin):
			return(tauMin);
		break;
		case(dbltauFactor):
			return(tauFactor);
		break;
		case(dblm):
			return(m);
		break;
		default:
			throw(std::invalid_argument("Unknown parameter"));
	}
}//end(Volume::get_dbl_par())
/*--------------------------------------------------------------------------*/

inline void Volume::SetVOLLog( std::ostream *outs , const char lvl ){
 Deflection::SetVOLLog( outs , lvl );

 #if( LOG_VOL )
  if( VOLLLvl > 1 )
   *VOLLog << std::endl <<  "Volume: tau = " << tauInit << " in ["<< tauMin
	   << " , -] with ( " << tauFactor << " , " << tauIter << " ) ~ m = "
	   << m << std::endl;
 #endif
 }
/*--------------------------------------------------------------------------*/

inline void Volume::NewDEF( void ){
 // change tau if needed - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
printf("Defl NrIter=%d\n", GetNrIter());
 if( GetNrIter() <= 1 ) {
  alpha = 1;
  return;
  }

 VarValue FiRec = ReadBestFiVal();//f(lambdabest)
printf("newdef bestfibar=%f\n", FiRec );
 // tau should be decreased if after tauIter iterations Fi() is decreased
 // less than 1% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! ( GetNrIter() % tauIter ) ){
	 printf("Dentro");
  if( tau > tauMin ) {
   VarValue x = ( lastvalue - FiRec ) / abs( FiRec );
   lastvalue = FiRec;
   if( x <= 0.01 ) {
    tau *= tauFactor;
    if( tau < tauMin )
     tau = tauMin;
    }
   }
 }
   printf("lastvalue=%f tau=%f\n", lastvalue, tau);
 // solve a quadratic problem and get the deflection coefficient - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 Solve();

 }  // end( Volume::NewDEF() )
/*--------------------------------------------------------------------------*/

inline c_VarValue Volume::GetDFLCoeff( void ){
 return( alpha );
 }
/*--------------------------------------------------------------------------*/

inline const VarValue Volume::Delta( void ){
	return( m * GetEpsLin() * std::max( 1.0 , std::abs( ReadBestFiVal() )) );
 }
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

inline int Volume::get_dflt_int_par_d( const idx_type par ){
	        	if( ( par >= inttauIter ) && (par<=intlast) )
	        	  return( dflt_int_par[par] );
	        	 else
	        	  throw(std::invalid_argument("Unknown parameter"));
 }
/*--------------------------------------------------------------------------*/

inline double Volume::get_dflt_dbl_par_d( const idx_type par ) {
	        	if( ( par >= dbltauInit ) && (par<=dblast) )
	        	  return( dflt_dbl_par[par] );
	        	 else
					throw(std::invalid_argument("Unknown parameter"));
}

/*--------------------------------------------------------------------------*/
/*-------------------------PRIVATE PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

inline void Volume::Solve( void ){
 VarValue GiNorm = GetGiNorm();
 VarValue DNorm = GetDNorm();
 VarValue GixD = GetdGk();
 VarValue Sigma = GetSigma();
 VarValue Epsilon = GetEpsilon();
 // compute the proximity parameter  - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
printf("GiNorm=%f, DNorm=%f, GixD=%f, Sigma=%f, Epsilon=%f\n", GiNorm, DNorm, GixD, Sigma, Epsilon);
 VarValue proximity = GetStepsize() == Inf<VarValue>() ? 1 : GetStepsize();

 // solving the cutting plane model  - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VarValue den = proximity * ( GiNorm + DNorm - 2.0 * GixD );

 if( den <= 1e-8 )
  alpha = 1;
 else
  alpha = ( Epsilon - Sigma - proximity * ( GixD - DNorm ) ) / den;

 // tau should be projected over [0,1] - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( alpha > 1 )
  alpha = std::min( tau , 1.0 );
 else
  if( alpha < 1e-8 )
   alpha = tau / 10;
printf("alpha=%f\n", alpha);
 }  // end ( void Solve )

/*--------------------------------------------------------------------------*/

}// end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Volume.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File Volume.h ------------------------------*/
/*--------------------------------------------------------------------------*/
