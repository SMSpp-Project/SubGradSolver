/*--------------------------------------------------------------------------*/
/*--------------------------- File Deflection.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the abstract class Deflection, which defines the interface
 * of the deflection rules (DR) used by the SubGradSolver class [see
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
 * \author Laura Fiorella \n
 *         Dipartimento di Matematica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Gorgone, Laura Fiorella
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __Deflection
 #define __Deflection /* self-identification: #endif at the end of the file */

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

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Deflection_CLASSES Classes in Deflection.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Deflection ------------------------------*/
/*--------------------------------------------------------------------------*/
/// interface of the deflection rules of SubGradSolver
/** The class Deflection provides the interface of the deflection rule used
 * by SubGradSolver [see SubGradSolver.h], whose task is to compute the
 * deflection coefficient \f$ \alpha_i \in ( 0 , 1 ] \f$ defining the
 * direction
 * \f[
 *   d_i = \alpha_i g_i + ( 1 - \alpha_i ) d_{i-1}
 * \f]
 * as a convex combination of the subgradient \f$ g_i \f$ and of the
 * previous direction. A rule may also decide whether the stability center
 * moves to the new point (a *serious step*) or not (a *null step*) [see
 * do_SS()]. No deflection at all, \f$ \alpha_i = 1 \f$, is obtained by
 * giving SubGradSolver no Deflection.
 *
 * The concrete rules register themselves in the Deflection factory with the
 * SMSpp_insert_in_factory_h / SMSpp_insert_in_factory_cpp_0 macros, and are
 * constructed by new_Deflection(). As for Stepsize, the parameters of the
 * rules belong to the SubGradSolver using them [see read_int_par() and
 * read_dbl_par()]. */

class Deflection
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

 /// type of the factory of Deflection
 using DeflectionFactory =
  boost::function< Deflection *( SubGradSolver * ) >;

 /// type of the map between strings and the factory of Deflection
 using DeflectionFactoryMap = std::map< std::string , DeflectionFactory >;

/*--------------------------------------------------------------------------*/
/*---------------------- CONSTRUCTOR AND DESTRUCTOR ------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor: takes the SubGradSolver using the rule
 explicit Deflection( SubGradSolver * slvr = nullptr ) : Solver( slvr ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual ~Deflection() = default;  ///< destructor: it is virtual, and empty

/*--------------------------------------------------------------------------*/
/*------------------------------ FACTORY -----------------------------------*/
/*--------------------------------------------------------------------------*/

 /// construct a Deflection of the given class for the given SubGradSolver
 /** Construct an object of the :Deflection class whose name is \p classname,
  * which must have been registered in the Deflection factory, for the
  * SubGradSolver \p slvr. Throws std::invalid_argument if no such class is
  * there. */

 static Deflection * new_Deflection( const std::string & classname ,
				     SubGradSolver * slvr ) {
  const std::string classname_( SMSpp_classname_normalise(
						 std::string( classname ) ) );
  const auto it = f_factory().find( classname_ );
  if( it == f_factory().end() )
   throw( std::invalid_argument( "Deflection::new_Deflection: " + classname
				 + " not present in Deflection factory" ) );
  return( ( it->second )( slvr ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the classname of the Deflection, i.e., its name in the factory

 [[nodiscard]] const std::string & classname( void ) const {
  return( private_name() );
  }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 /// (re)initialize the rule; it does nothing in the base class
 virtual void format( void ) {}

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR DEFLECTION COMPUTATION ------------------*/
/*--------------------------------------------------------------------------*/

 /// compute the deflection coefficient of the current iteration
 /** This method must be called before get_deflection(), and it is the core
  * of every derived class. */

 virtual void new_deflection( void ) = 0;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// returns the deflection coefficient computed by new_deflection()
 [[nodiscard]] virtual FunctionValue get_deflection( void ) const = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// tells whether a serious step takes place
 /** Returns true if the stability center has to move to the new point (a
  * *serious step*, SS), false if it stays where it is (a *null step*, NS).
  * Typically, a SS occurs when the decrease of the function is as good as
  * the expected one [see delta()]. The default is true. */

 [[nodiscard]] virtual bool do_SS( void ) { return( true ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the expected decrease of the function [see do_SS()]
 /** Returns the expected decrease of the function, if the rule has such a
  * concept; the default is + INF. */

 [[nodiscard]] virtual FunctionValue delta( void ) const {
  return( Inf< FunctionValue >() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// tells whether the rule also chooses the stepsize
 /** A rule that chooses both the stepsize and the deflection coefficient,
  * as PrimalDual does, can only be given as the Stepsize [see
  * Stepsize::deflection_role()]; the default is false. */

 [[nodiscard]] virtual bool is_also_stepsize( void ) const {
  return( false );
  }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/** @name Protected methods to read from SubGradSolver
 *
 * Deflection is a friend of SubGradSolver, but the classes derived from it
 * are not: these methods give them what they need of SubGradSolver. They are
 * implemented in SubGradSolver.cpp.
 *  @{ */

 /// the stepsize \f$ \nu_i \f$ of the last iteration
 [[nodiscard]] FunctionValue read_stepsize( void ) const;

 /// the squared norm of the subgradient \f$ g_i \f$
 [[nodiscard]] FunctionValue read_Gi_norm( void ) const;

 /// the squared norm of the direction \f$ d_{i-1} \f$
 [[nodiscard]] FunctionValue read_dir_norm( void ) const;

 /// the scalar product \f$ d_{i-1}^\top g_i \f$
 [[nodiscard]] FunctionValue read_dGk( void ) const;

 /// the linearization error \f$ \sigma_i \f$ of \f$ g_i \f$ at the center
 [[nodiscard]] FunctionValue read_sigma( void ) const;

 /// the linearization error \f$ \epsilon_{i-1} \f$ of \f$ d_{i-1} \f$
 [[nodiscard]] FunctionValue read_epsilon( void ) const;

 /// the value at the last evaluated point \f$ \lambda_i \f$
 [[nodiscard]] FunctionValue read_Fi_lambda( void ) const;

 /// the value at the stability center \f$ \bar{\lambda}_i \f$
 [[nodiscard]] FunctionValue read_Fi_bar( void ) const;

 /// the best value found so far
 [[nodiscard]] FunctionValue read_Fi_best( void ) const;

 /// the total number of iterations since the rule has been formatted
 [[nodiscard]] Index read_iterations( void ) const;

 /// the (estimate of the) Lipschitz constant of the function
 [[nodiscard]] FunctionValue read_Lipschitz( void ) const;

 /// the parameter dbltStar of SubGradSolver
 [[nodiscard]] FunctionValue read_tStar( void ) const;

 /// the relative accuracy dblRelAcc of SubGradSolver
 [[nodiscard]] FunctionValue read_rel_acc( void ) const;

 /// an int parameter of SubGradSolver
 [[nodiscard]] int read_int_par( idx_type par ) const;

 /// a double parameter of SubGradSolver
 [[nodiscard]] double read_dbl_par( idx_type par ) const;

/** @} ---------------------------------------------------------------------*/
/*---------------------------- FACTORY -------------------------------------*/
/*--------------------------------------------------------------------------*/

 /// the factory of Deflection
 static DeflectionFactoryMap & f_factory( void );

 /// empty placeholder for class-specific static initialization
 static void static_initialization( void ) {}

 /// the name of the class, as registered in the factory
 [[nodiscard]] virtual const std::string & private_name( void ) const = 0;

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

 SubGradSolver * Solver;  ///< (pointer to) the SubGradSolver

/*--------------------------------------------------------------------------*/

 };  // end( class Deflection )

/** @} end( group( Deflection_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Deflection.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File Deflection.h ---------------------------*/
/*--------------------------------------------------------------------------*/
