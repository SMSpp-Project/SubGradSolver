/*--------------------------------------------------------------------------*/
/*--------------------------- File Stepsize.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the abstract class Stepsize, which defines the interface of
 * the stepsize rules (SR) used by the SubGradSolver class [see
 * SubGradSolver.h], together with the factory the concrete rules register
 * themselves in, so that a SubGradSolver can select its rule by name.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Matematica ed Informatica \n
 *         Universita' di Cagliari \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Gorgone, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __Stepsize
 #define __Stepsize /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "C05Function.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

class SubGradSolver;     // forward declaration of class SubGradSolver
class Deflection;        // forward declaration of class Deflection

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Stepsize_CLASSES Classes in Stepsize.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Stepsize --------------------------------*/
/*--------------------------------------------------------------------------*/
/// interface of the stepsize rules of SubGradSolver
/** The class Stepsize provides the interface of the stepsize rule used by
 * SubGradSolver [see SubGradSolver.h]. The main methods are new_step(),
 * which sets up the stepsize of the current iteration, and get_stepsize(),
 * which returns it.
 *
 * The base class implements the formula of the *target value* stepsize
 * rules, i.e.,
 * \f[
 *   \nu_i = \beta_i \frac{ f( \bar{\lambda}_i ) - f^{lev}_i }{ \| d_i \|^2 }
 * \f]
 * where \f$ \bar{\lambda}_i \f$ is the stability center, \f$ d_i \f$ the
 * direction, \f$ \beta_i \in ( 0 , 2 ) \f$ and \f$ f^{lev}_i \f$ the target
 * level; the derived classes then only have to define how \f$ \beta_i \f$
 * and \f$ f^{lev}_i \f$ are chosen in new_step(). All the values refer to the
 * function that SubGradSolver minimizes, i.e., to the opposite of the
 * objective when this is a concave one to be maximized.
 *
 * The concrete rules register themselves in the Stepsize factory with the
 * SMSpp_insert_in_factory_h / SMSpp_insert_in_factory_cpp_0 macros, and are
 * constructed by new_Stepsize(). The parameters of the rules belong to the
 * SubGradSolver using them, which is what allows to set them in the very
 * same ComputeConfig of the SubGradSolver, even before the rule is
 * constructed: a rule reads them from its SubGradSolver when it needs them
 * [see read_int_par() and read_dbl_par()]. */

class Stepsize
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 using Index = Function::Index;                     ///< an index
 using FunctionValue = Function::FunctionValue;     ///< a function value
 using idx_type = ThinComputeInterface::idx_type;   ///< a parameter index

 /// type of the factory of Stepsize
 using StepsizeFactory = boost::function< Stepsize *( SubGradSolver * ) >;

 /// type of the map between strings and the factory of Stepsize
 using StepsizeFactoryMap = std::map< std::string , StepsizeFactory >;

/*--------------------------------------------------------------------------*/
/*---------------------- CONSTRUCTOR AND DESTRUCTOR ------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor: takes the SubGradSolver using the rule
 explicit Stepsize( SubGradSolver * slvr = nullptr ) : Solver( slvr ) {
  FiLev = - Inf< FunctionValue >();
  Beta = MaxBeta = Inf< FunctionValue >();
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual ~Stepsize() = default;  ///< destructor: it is virtual, and empty

/*--------------------------------------------------------------------------*/
/*------------------------------ FACTORY -----------------------------------*/
/*--------------------------------------------------------------------------*/

 /// construct a Stepsize of the given class for the given SubGradSolver
 /** Construct an object of the :Stepsize class whose name is \p classname,
  * which must have been registered in the Stepsize factory, for the
  * SubGradSolver \p slvr. Throws std::invalid_argument if no such class is
  * there. */

 static Stepsize * new_Stepsize( const std::string & classname ,
				 SubGradSolver * slvr ) {
  const std::string classname_( SMSpp_classname_normalise(
						 std::string( classname ) ) );
  const auto it = f_factory().find( classname_ );
  if( it == f_factory().end() )
   throw( std::invalid_argument( "Stepsize::new_Stepsize: " + classname +
				 " not present in Stepsize factory" ) );
  return( ( it->second )( slvr ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the classname of the Stepsize, i.e., its name in the factory

 [[nodiscard]] const std::string & classname( void ) const {
  return( private_name() );
  }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 /// (re)initialize the rule
 /** Initializes the rule, forgetting everything it has learnt about the
  * function: it is called by SubGradSolver before the first iteration, and
  * every time the function changes in a way that makes that information
  * useless. */

 virtual void format( void ) {
  Beta = MaxBeta = Inf< FunctionValue >();
  FiLev = - Inf< FunctionValue >();
  }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
/*--------------------------------------------------------------------------*/

 /// set up the stepsize of the current iteration
 /** This method must be called before get_stepsize(), and it is the core of
  * every derived class. A target value stepsize rule computes here both
  * \f$ \beta_i \f$ and \f$ f^{lev}_i \f$. */

 virtual void new_step( void ) = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// tells whether the rule needs the scalar product of d_{i-1} and g_i
 /** Returns true if the rule needs the scalar product \f$ d_{i-1}^\top g_i
  * \f$, which SubGradSolver then computes; the default is false. */

 [[nodiscard]] virtual bool needs_dM1Gk( void ) const { return( false ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the Deflection role of the rule, if it also has one
 /** A rule that chooses both the stepsize and the deflection coefficient,
  * as PrimalDual does, returns itself as a Deflection, which SubGradSolver
  * then uses in that role too; the default is nullptr. */

 [[nodiscard]] virtual Deflection * deflection_role( void ) {
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// returns the stepsize of the current iteration
 /** Returns the stepsize computed by the last call to new_step(). The base
  * class implements the target value formula; for the incremental steps
  * (\p incremental == true) of the iterations \f$ k = p( \Pi + 1 ) , \ldots
  * , ( p + 1 )( \Pi + 1 ) - 1 \f$ the formula is
  * \f[
  *   \nu_k = \beta_k \frac{ f_{ p( \Pi + 1 ) } - f^{lev}_{ p( \Pi + 1 ) } }
  *                        { \chi ( \Pi + 1 ) L^2 }
  * \f]
  * where \f$ \chi \f$ is dblLipFactor, \f$ L \f$ the Lipschitz constant of
  * the function and \f$ \Pi \f$ the number of incremental iterations. When
  * the denominator is (almost) zero the step is 1. */

 virtual FunctionValue get_stepsize( bool incremental = false );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the target level of the current iteration

 [[nodiscard]] virtual FunctionValue get_level( void ) const {
  return( FiLev );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the factor beta of the current iteration

 [[nodiscard]] virtual FunctionValue get_beta( void ) const {
  return( Beta );
  }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

 /// change the maximum value of beta
 /** Changes the maximum value of \f$ \beta_i \f$ of a target value rule,
  * which is how SubGradSolver implements the safe rule of the
  * stepsize-restricted scheme, by passing the deflection coefficient
  * \f$ \alpha_i \f$. */

 virtual void set_max_beta( FunctionValue alpha ) { MaxBeta = alpha; }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// tries to improve the target level, returns true if it has changed
 virtual bool update_target_level( void ) { return( false ); }

/*--------------------------------------------------------------------------*/
/** @name Protected methods to read from SubGradSolver
 *
 * Stepsize is a friend of SubGradSolver, but the classes derived from it are
 * not: these methods give them what they need of SubGradSolver. They are
 * implemented in SubGradSolver.cpp.
 *  @{ */

 /// the deflection coefficient \f$ \alpha_i \f$
 [[nodiscard]] FunctionValue read_alpha( void ) const;

 /// the squared norm of the subgradient \f$ g_i \f$
 [[nodiscard]] FunctionValue read_Gi_norm( void ) const;

 /// the squared norm of the direction \f$ d_i \f$
 [[nodiscard]] FunctionValue read_dir_norm( void ) const;

 /// the scalar product \f$ d_i^\top g_i \f$
 [[nodiscard]] FunctionValue read_dGk( void ) const;

 /// the scalar product \f$ d_{i-1}^\top g_i \f$
 [[nodiscard]] FunctionValue read_dM1Gk( void ) const;

 /// the number of incremental iterations after an outer one
 [[nodiscard]] Index read_NItIncr( void ) const;

 /// the value at the last evaluated point \f$ \lambda_i \f$
 [[nodiscard]] FunctionValue read_Fi_lambda( void ) const;

 /// the value at the stability center \f$ \bar{\lambda}_i \f$
 [[nodiscard]] FunctionValue read_Fi_bar( void ) const;

 /// the best value found so far
 [[nodiscard]] FunctionValue read_Fi_best( void ) const;

 /// the best known (possibly conditional) lower bound, - INF if none
 [[nodiscard]] FunctionValue read_lower_bound( void ) const;

 /// the (estimate of the) Lipschitz constant of the function
 [[nodiscard]] FunctionValue read_Lipschitz( void ) const;

 /// the parameter dbltStar of SubGradSolver
 [[nodiscard]] FunctionValue read_tStar( void ) const;

 /// an int parameter of SubGradSolver
 [[nodiscard]] int read_int_par( idx_type par ) const;

 /// a double parameter of SubGradSolver
 [[nodiscard]] double read_dbl_par( idx_type par ) const;

/** @} ---------------------------------------------------------------------*/
/*---------------------------- FACTORY -------------------------------------*/
/*--------------------------------------------------------------------------*/

 /// the factory of Stepsize
 static StepsizeFactoryMap & f_factory( void );

 /// empty placeholder for class-specific static initialization
 static void static_initialization( void ) {}

 /// the name of the class, as registered in the factory
 [[nodiscard]] virtual const std::string & private_name( void ) const = 0;

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

 SubGradSolver * Solver;  ///< (pointer to) the SubGradSolver

 FunctionValue FiLev;     ///< the target level \f$ f^{lev}_i \f$
 FunctionValue Beta;      ///< the factor \f$ \beta_i \f$
 FunctionValue MaxBeta;   ///< the maximum value for \f$ \beta_i \f$

/*--------------------------------------------------------------------------*/

 };  // end( class Stepsize )

/** @} end( group( Stepsize_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Stepsize.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File Stepsize.h -----------------------------*/
/*--------------------------------------------------------------------------*/
