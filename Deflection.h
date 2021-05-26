/*--------------------------------------------------------------------------*/
/*--------------------------- File Deflection.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Description of the abstract class <tt>Deflection</tt>, which defines the
 * interface for computing the delflection coefficient within the SubGradSolver
 * class [see SubGradSolver.h].
 * 
 * \version 0.1
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

#ifndef __Deflection
 #define __Deflection /* self-identification: #endif at the end of the file */


/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/
namespace SMSpp_di_unipi_it
{

 class SubGradSolver;     // forward declaration of class SubGradSolver

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Deflection ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Deflection_CLASSES Classes in Deflection.h
    @{ */

/** The class Deflection provides an interface for the deflection rule (DR),
 *  to be used in the SubGradSolver solver [see SubGradSolver.h].
 *
 *  The aim of this class is to compute the deflection coefficient
 *  \f$ \alpha_i \f$, which is employed in the search direction formula as
 *  follows:
 *  \f[
 *   d_i = \alpha_i g_i + (1-\alpha_i)d_{i-1}\;.
 *  \f]
 *
 *  The user must extend the class to one or more deflection rules following
 *  this interface. */

class Deflection
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
 
 using VarValue = Function::FunctionValue;
 using c_VarValue = Function::c_FunctionValue;

 using Vec_VarValue = Function::Vec_FunctionValue;
using Index = Function::Index;
using idx_type=ThinComputeInterface::idx_type;
/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{ */

/** Constructor of the class. It has no parameters apart the SubGradSolver object
 * using it because there is no "standard" deflection rule (except "no
 * deflection", i.e., \f$ \alpha_i = 1 \f$, which however in SubGradSolver is
 * treated by just setting the null Deflection object). */

   Deflection( SubGradSolver * slvr )
   {
    if( ! slvr )
    throw (std::invalid_argument( "Deflection: no subgradient solver " ));
    Solver = slvr;

    VOLLog = 0;
    VOLLLvl = 0;

    }  // end( Deflection::Deflection )

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

/** The class outputs "log" information onto the ostream pointed by outs.
 * lvl controls the "level of verbosity" of the code: lvl == 0 means that
 * nothing at all is printed, and values larger than 0 mean increasing
 * amounts of information, the specific effect of each value being derived-
 * class-dependent. outs == 0 implies lvl == 0. */

   virtual void SetVOLLog( std::ostream *outs = 0 , const char lvl = 0 )
   {
    if( ( VOLLog = outs ) )
     VOLLLvl = lvl;
    else
     VOLLLvl = 0;
    }

/*--------------------------------------------------------------------------*/
/** The method initializes the DR. It does nothing for the base class. */

   virtual void Format( void ) {};

/**@} ----------------------------------------------------------------------*/
/*-------------------- METHODS FOR DEFLECTION COMPUTATION ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Deflection computation
    @{ */

/** This method must be called before GetDFLCoeff() [see below].
 * This is indeed the core of every derived class, producing a new
 * deflection coefficient.
 *
 * Typically, the previous coefficient will be unavailable after the call to
 * NewDEF() [GetDFLCoeff()]. */

   virtual void NewDEF( void ) = 0;

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
    @{ */

/** Returns the deflection coefficient. This function must be called after
   NewDEF() [see above]. */

   virtual c_VarValue GetDFLCoeff( void ) = 0;

/*--------------------------------------------------------------------------*/
/** Returns true if a <em>Serious Step</em> (SS) takes place. Otherwise a
 * <em>Null Step (NS)</em> occurs. By default, true is returned.
 *
 * Typically, a SS occurs when a \" good \" improvement of the function
 * \f$ f \f$ is obtained, usually it should be as good as the expected 
 * improvement [see Delta()]. */

   virtual const bool DoSS( void ) { return( true ); }

/*--------------------------------------------------------------------------*/
/** Returns the expected improvement in the objective function [see DoSS()].
 * This value is not useful for all DR. By default Inf<VarValue>() is returned.
 */

   virtual c_VarValue Delta( void ) { return( Inf<VarValue>() ); }

/**@} ----------------------------------------------------------------------*/
/*------------------------------ ACCESS TO PARAMETERS----------------------*/

virtual void set_par_d(const idx_type par , const int value ){};

virtual void set_par_d(const idx_type par , const double value){};

virtual int get_int_par_d(const idx_type par)=0;

virtual double get_dbl_par_d (const idx_type par)=0;

virtual int get_dflt_int_par_d( const idx_type par )=0;

virtual double get_dflt_dbl_par_d( const idx_type par )=0;

/*--------------------------------------------------------------------------*/

  static Deflection * new_Deflection( const std::string & classname, SubGradSolver * slvr = nullptr ) {
  const std::string classname_( SMSpp_classname_normalise(
					        std::string( classname ) ) );
  const auto it = Deflection::f_factory().find( classname_ );
  if( it == Deflection::f_factory().end() )
   throw( std::invalid_argument( classname +
				 " not present in Deflection factory" ) );
  return( ( it->second )( slvr ) );
  }
  
   const std::string & classname( void ) { return( private_name() ); }
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~Deflection() {}

/**@} ----------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:
 
 /*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods to read from SubGradSolver
 *
 * These methods are used to provide access to the information stored into
 * the SubGradSolver object. Deflection is a "friend" of SubGradSolver and 
 * therefore it can read its protected data structures, but classes derived 
 * from Deflection are not friend of SubGrad and cannot. This is why these 
 * methods are defined in the base class (and implemented in SubGradSolver.cpp).
 * @{ */

/*--------------------------------------------------------------------------*/
/** Returns the stepsize \f$\nu_i\f$. */

   VarValue GetStepsize( void );

/*--------------------------------------------------------------------------*/
/** Returns the norm of the subgradient \f$ g_i \f$. */

   VarValue GetGiNorm( void );

/*--------------------------------------------------------------------------*/
/** Returns the norm of the direction \f$ d_i \f$. */

   VarValue GetDNorm( void );

/*--------------------------------------------------------------------------*/
/** Returns the scalar product \f$ g_i^{\top} d_i \f$. */

   VarValue GetdGk( void );

/*--------------------------------------------------------------------------*/
/** Returns the linearization error \f$ \sigma_i \f$. */

   VarValue GetSigma( void );

/*--------------------------------------------------------------------------*/
/** Returns the linearization error \f$ \epsilon_i \f$. */

   VarValue GetEpsilon( void );

/*--------------------------------------------------------------------------*/
/** Returns the full function \f$ f(\lambda_i) \f$. */

   VarValue ReadFVal( void );
   
/*--------------------------------------------------------------------------*/
/** Returns the full function \f$ f(\bar\lambda_i) \f$. */

   VarValue ReadFiVal( void );

/*--------------------------------------------------------------------------*/
/** Returns the number of iteration done this far. */

   Index GetNrIter( void );

/*--------------------------------------------------------------------------*/
/** Returns the the best value of \f$ f \f$ found this far. */

   VarValue ReadBestFiVal( void );

/*--------------------------------------------------------------------------*/
/** Returns the Lipschitz constant of the objective function. */

   VarValue GetGlobalLipschitz( void );

/*--------------------------------------------------------------------------*/
/** Returns the default  error. */

   VarValue GetEpsLin( void );

/*--------------------------------------------------------------------------*/
/**@}-----------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

  SubGradSolver *Solver;    ///< (pointer to) the SubGradSolver object
  std::ostream *VOLLog;    ///< the output stream object
  char VOLLLvl;       ///< the "level of verbosity"
  
 
/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/
/**
 * Factory definition
 **/
typedef boost::function< Deflection * (SubGradSolver *) > DeflectionFactory;
// type of the factory of Deflection

typedef std::map< std::string, DeflectionFactory > DeflectionFactoryMap;
// type of the map between strings and the factory of Deflection

static DeflectionFactoryMap & f_factory( void );

static void static_initialization(void){};
 
/*--------------------------------------------------------------------------*/
private:

virtual const std::string & private_name( void ) const = 0;

 };  // end( class Deflection )

/** @} end( group( Deflection_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

};  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif 

/*--------------------------------------------------------------------------*/
/*------------------------ End File Deflection.h ---------------------------*/
/*--------------------------------------------------------------------------*/
