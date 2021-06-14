/*--------------------------------------------------------------------------*/
/*--------------------------- File Stepsize.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Definition of the abstract <tt>Stepsize</tt> class, which defines the
 * interface for computing the stepsize rule (SR) used by the
 * <tt>SubGradSolver</tt> class [see SubGradSolver.h].
 *
 * \version 1.02
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __Stepsize
 #define __Stepsize /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{

class SubGradSolver;     // forward declaration of class SubGradSolver

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Stepsize --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Stepsize_CLASSES Classes in Stepsize.h
    @{ */

/** The class Stepsize provides an interface for the stepsize rule, to be
 *  used in the SubGradSolver [see SubGradSolver.h].
 *
 *  The main functions of the class are NewStep() and GetStepsize(). While
 *  NewStep() sets up a new step, the value \f$\nu_i\f$ is passed to the
 *  SubGradSolver by means of a call to GetStepsize().
 *
 *  The user must extend the class to one or more stepsize rules following
 *  this interface.  */

class Stepsize
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
 /*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC TYPES OF THE CLASS ------------------------*/
enum step_dbl_par{
	dblLpsFct=0,
	lastdblstep,
};

// define and initialize here the vector of double parameters names
const std::vector< std::string > dbl_par_string={
	"dblLpsFct"};

// define and initialize here the map for double parameters names
const std::map< std::string , ThinComputeInterface::idx_type > dbl_pars_map ={
	{"dblLpsFct", dblLpsFct}};

// define and initialize here the default double parameters
const std::vector<double> dflt_dbl_par ={
	1
};

 using VarValue = Function::FunctionValue;
 using c_VarValue = Function::c_FunctionValue;

 using Vec_VarValue = Function::Vec_FunctionValue;
 
 using Index = Function::Index;
 using c_Index = Function::c_Index;

 using idx_type =ThinComputeInterface::idx_type;
/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{  Only a few parameters are deemed to be general enough to be put in the
 * basic Stepsize class:
 * 		VarValue LpsFct [1] a general scaling factor (for instance, it is the 
 *                  		Lipschitz factor in the target value SR)  */

   Stepsize( SubGradSolver *slvr)
   {
    if( ! slvr )
     throw (std::invalid_argument("Stepsize: no subgradient solver "));
		Solver = slvr;
    LpsFct=dflt_dbl_par[dblLpsFct];
    STPLog = 0;
    STPLLvl = 0;

    }  // end( Stepsize::Stepsize() )

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

/** The class outputs "log" information onto the ostream pointed by outs.
   lvl controls the "level of verbosity" of the code: lvl == 0 means that
   nothing at all is printed, and values larger than 0 mean increasing
   amounts of information, the specific effect of each value being derived-
   class-dependent. outs == 0 implies lvl == 0. */

   virtual void SetSTPLog( std::ostream *outs = 0 , const char lvl = 0 )
   {
    if( ( STPLog = outs ) )
     STPLLvl = lvl;
    else
     STPLLvl = 0;
    }

/*--------------------------------------------------------------------------*/
/** The method initializes the stepsize rule. */

   virtual void Format( void )
   {
    Beta = MaxBeta = Inf<VarValue>();
    FiLev = - Inf<VarValue>();
    }

/**@} ----------------------------------------------------------------------*/
/*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Stepsize computation
    @{ */

/** This method must be called before GetStepsize() [see below]. This is
   indeed the core of every derived class, producing a new step.

   If the derived class implements a target value stepsize rule, the method
   must compute both the scalar \f$ \beta_i \f$ and the level
   \f$ f^{lev}_i \f$.

   Typically, the previous step will be unavailable after the call to
   NewStep() [see GetStepsize()].  */

   virtual void NewStep( void ) = 0;

/*--------------------------------------------------------------------------*/
/** Returns true if the stepsize requires the scalar product
   \f$ d_{i-1}^{\top} g_i\f$. The default implementation don't require this
   value. */

   virtual bool NeedsdkM1Gk( void ) { return( false ); }

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
    @{ */

/** This function must be called after NewStep() [see above ] and returns the
 * step value. StepIsIncr controls the stepsize formula in the incremental
 * version.
 *
 * The base class provides the step value for target stepsize rules, by
 * exploiting \f$ beta_k \f$ and \f$ f^{lev}_k \f$ computed by NewStep().
 * While for non-incremental step (StepIsIncr = false) for \f$ k = 1, 2 ,
 * \ldots \f$ the formula is
 * \f[
 *   \nu_k = \beta_k * ( f_k - f^{lev}_k ) / \| d_k \|^2\,,
 * \f]
 * for incremental step (StepIsIncr = true) for \f$ k = p(\Pi+1) , \ldots ,
 * (p+1)(\Pi+1) - 1 , p = 0, 1 , 2 \ldots \f$ the stepsize rule is
 * \f[
 *   \nu_k = \beta_k \frac{ f_{ p(\Pi + 1 ) } - f^{lev}_{ p( \Pi + 1 ) } }
 *                        { \chi( \Pi + 1 ) C^2 }
 * \f]
 * where \f$ \chi \f$ is the LpsFct factor and \f$ C \f$ is the Lipschitz
 * constant.
 *
 * Note that GetStepsize() only returns the value, the actual step
 * computation being done by NewStep(). */
 

   virtual VarValue GetStepsize( bool StepIsIncr = false ){
    VarValue denomin, step;
    const VarValue STEpZro = 1e-8;
    printf("StepIsIncr=%d\n", StepIsIncr);
    if( StepIsIncr ) {
     VarValue LipCns = GetGlobalLipschitz();
     printf("LipCns=%f\n", LipCns);
     denomin = LipCns * LipCns * VarValue( LpsFct * GetNItIcr() );
     }
    else
     denomin = GetDNorm();
    printf("denomin=%f e FiLev=%f e STEpZro=%f e FiBar=%f e Beta=%f e Maxbeta=%f\n", denomin, FiLev, STEpZro, ReadFiBar(), Beta, MaxBeta);
    if( denomin <= ( STEpZro * std::max( abs( ReadFiBar() ) , VarValue(1) ) ) )
     step = 1;
    else {
     VarValue Beta_ = Beta;
     if( Beta_ >= std::min( 2.0 , MaxBeta ) )
      Beta_ = std::min( 2.0 - STEpZro , MaxBeta );
      printf("Beta_=%f\n", Beta_);
     step = Beta_ * ( ReadFiBar() - FiLev ) / denomin;
     }
     printf("denomin=%f e FiLev=%f e STEpZro=%f e FiBar=%f e Beta=%f e Maxbeta=%f\n", denomin, FiLev, STEpZro, ReadFiBar(), Beta, MaxBeta);
    printf("step=%f\n", step);
    return( step );

    }  // end( GetStepsize() )

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR ALGORITHMIC PARAMETERS -------------------*/
/*--------------------------------------------------------------------------*/

  virtual void set_par(const idx_type par , const int value );

  virtual void set_par(const idx_type par , const double value);

  virtual int get_int_par(const idx_type par) ;

  virtual double get_dbl_par(const idx_type par) ;

  virtual double get_dflt_dbl_par( const idx_type par ) ;

  virtual int get_dflt_int_par( const idx_type par ) ;
    
/*--------------------------------------------------------------------------*/
/** This function must be called after NewStep() [see above].
 *
 * The default implementation returns the level \f$ f^{lev}_i \f$ for a
 * target value stepsize rule [see below]. It should throw exception for
 * other forms of the stepsize rule. */

   virtual VarValue GetLev( void ) { return( FiLev ); }

/*--------------------------------------------------------------------------*/
/** This function must be called after NewStep() [see above].
 *
 * The default implementation returns the scalar \f$ \beta_i \f$ for a
 * target value stepsize rule [see below]. It should throw exception for
 * other forms of the stepsize rule. */

   virtual VarValue GetBeta( void ) { return( Beta ); }

/**@} ----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
    @{ */

/** Changes the <em>maximum value</em> of \f$ \beta_i \f$ for a <em>target
 * value stepsize rule</em>.
 *
 * Typically, the method implements the safe rule of the stepsize-restricted
 * approach. Passing the deflection coefficient \f$ \alpha_i \f$ does exactly
 * this job. */

   virtual void SetMaxBeta( const VarValue alpha )
   {
    MaxBeta = alpha;
    }

/**@} ----------------------------------------------------------------------*/

static Stepsize * new_Stepsize( const std::string & classname,
 SubGradSolver * slvr = nullptr ) {
  const std::string classname_( SMSpp_classname_normalise(
					        std::string( classname ) ) );
  const auto it = Stepsize::f_factory().find( classname_ );
  if( it == Stepsize::f_factory().end() )
   throw( std::invalid_argument( classname +
				 " not present in Stepsize factory" ) );
  return( ( it->second )(slvr) );
  }


 const std::string & classname( void ) { return( private_name() ); }
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~Stepsize() {}

/**@} ----------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

  protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** It tries to provide a new target level \f$ f^{lev}_i \f$. It true is
 * returned, the level has been changed. By default, it returns false. */

   virtual bool UpdateTargetLevel( void ) { return( false ); }

/*--------------------------------------------------------------------------*/
/** @name Protected methods to read from SubgradSolver
 *
 * These methods are used to provide access to the information stored into
 * the SubGradSolver object. Stepsize is a "friend" of SubGradSolver and 
 * therefore it can read its protected data structures, but classes derived 
 * from Stepsize are not friend of SubGradSolver and cannot. 
 * This is why these methods are defined in the base class 
 * (and implemented in SubGradSolver.cpp).
 * @{ */

/*--------------------------------------------------------------------------*/
/** Returns the deflection coefficient \f$\alpha_i\f$. */

  VarValue GetCoeffDefl( void );

/*--------------------------------------------------------------------------*/
/** Returns the norm of the subgradient \f$ g_i \f$. */

  VarValue GetGiNorm( void );

/*--------------------------------------------------------------------------*/
/** Returns the norm of the direction \f$ d_i \f$. */

  VarValue GetDNorm( void );

/*--------------------------------------------------------------------------*/
/** Returns the scalar product \f$ g_i^{\top} d_i\f$. */

  VarValue GetdGk( void );

/*--------------------------------------------------------------------------*/
/** Returns the scalar product \f$ g_i^{\top} d_{i-1}\f$. */

  VarValue GetdkM1Gk( void );

/*--------------------------------------------------------------------------*/
/** Returns NItIncr, the parameter for incremental iterations [see SubGrad.h].
 */
  Index GetNItIcr( void );

/*--------------------------------------------------------------------------*/
/** Returns FiLambda [ see SubGrad.h ]. */

   VarValue ReadFkVal( void );
    
/*--------------------------------------------------------------------------*/
/** Returns FiBar. */
    
   VarValue ReadFiBar( void );

/*--------------------------------------------------------------------------*/
/** Returns the the best value of \f$ f \f$ found this far. */

   VarValue ReadBestFiVal( void );

/*--------------------------------------------------------------------------*/
/** Returns the Lipschitz constant of the objective function. */

   VarValue GetGlobalLipschitz( void );

/*--------------------------------------------------------------------------*/
/** Returns the Lowerbound. */

   VarValue GetLowerbound( void );

/*--------------------------------------------------------------------------*/
/** Returns the longest step can be performed. */

   VarValue GettStar( void );

/*--------------------------------------------------------------------------*/

/**@} ----------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/
  
  SubGradSolver *Solver;  ///< (pointer to) the SubGradSolver solver
  std::ostream *STPLog;  ///< the output stream object
  char STPLLvl;     ///< the "level of verbosity"

  VarValue FiLev;      ///< the target level \f$ f^{lev}_i \f$
  VarValue Beta;       ///< beta factor \f$ \beta_i \f$
  VarValue MaxBeta;    ///< maximum value for beta factor

  VarValue LpsFct;     ///< scaling factor

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

typedef boost::function< Stepsize *(SubGradSolver *) > StepsizeFactory;
// type of the factory of Stepsize

typedef std::map< std::string, StepsizeFactory > StepsizeFactoryMap;
// type of the map between strings and the factory of Stepsize

static StepsizeFactoryMap & f_factory( void );

static void static_initialization(void){};



/*------------------------PRIVATE PART OF THE CLASS -----------------------*/
private:

virtual const std::string & private_name( void ) const = 0;

 }; // end( class Stepsize )

/** @} end( group( Stepsize_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

};  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Stepsize.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File Stepsize.h -----------------------------*/
/*--------------------------------------------------------------------------*/
