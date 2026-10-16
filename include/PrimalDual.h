/*--------------------------------------------------------------------------*/
/*--------------------------- File PrimalDual.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the class PrimalDual, which implements two variants of the
 * Primal-Dual subgradient method (PDSM), with simple and weighted averages.
 * These boil down to the *simultaneous* choice of a deflection coefficient
 * and of a stepsize, hence the class implements *both* the Deflection [see
 * Deflection.h] and the Stepsize [see Stepsize.h] interfaces. For details
 * see
 *
 *  Y. Nesterov. "Primal-dual subgradient methods for convex problems"
 *  Mathematical Programming 120, 221-259, 2009
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

#ifndef __PrimalDual
 #define __PrimalDual /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Deflection.h"
#include "Stepsize.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup PrimalDual_CLASSES Classes in PrimalDual.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS PrimalDual ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Primal-Dual subgradient method, both a Deflection and a Stepsize
/** The class PrimalDual implements the simple (intPDAverage == 1) and the
 * weighted (intPDAverage == 0) averages variants of the Primal-Dual
 * subgradient method in the unified scheme of SubGradSolver. The direction
 * is the average of all the previous weighted subgradients,
 * \f[
 *   d_i = \frac{ \sum_{k=1}^i \upsilon_k g_k }{ \Delta_i } \quad , \quad
 *   \Delta_i = \sum_{k=1}^i \upsilon_k
 * \f]
 * where \f$ \upsilon_k = 1 \f$ for simple averages and \f$ \upsilon_k = 1 /
 * \| g_k \| \f$ for weighted ones; the deflection coefficient and the
 * stepsize then are
 * \f[
 *   \alpha_i = \upsilon_i / \Delta_i \quad , \quad
 *   \nu_i = \Delta_i / \omega_i
 * \f]
 * where \f$ \omega_i = \gamma \hat{\beta}_i \f$, \f$ \hat{\beta}_1 = 1 \f$,
 * \f$ \hat{\beta}_{i+1} = \hat{\beta}_i + 1 / \hat{\beta}_i \f$, and
 * \f$ \gamma \f$ depends on the Lipschitz constant \f$ L \f$ of the
 * function, on dbltStar and on the factor F = dblLipFactor of SubGradSolver
 * as \f$ \gamma = F \sqrt{ L / 2 } / t^* \f$ for simple averages and
 * \f$ \gamma = F / ( \sqrt{ 2 L } t^* ) \f$ for weighted ones.
 *
 * The method takes all its steps from the first point, which is the prox
 * center of Nesterov's dual averaging, hence the stability center never
 * moves [see do_SS()].
 *
 * Being both a Deflection and a Stepsize, PrimalDual is registered in both
 * factories, and SubGradSolver uses one single object for both roles. Since
 * the deflection coefficient determines the stepsize, the deflection is
 * computed first. */

class PrimalDual : public Deflection , public Stepsize
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 using FunctionValue = Stepsize::FunctionValue;  ///< a function value

/*--------------------------------------------------------------------------*/
/*---------------------- CONSTRUCTOR AND DESTRUCTOR ------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor: takes the SubGradSolver using the rule
 explicit PrimalDual( SubGradSolver * slvr = nullptr )
  : Deflection( slvr ) , Stepsize( slvr ) , delta_( 0 ) , tmp_delta( 0 ) ,
    vk( 1 ) , step_( 0 ) , gamma( 1 ) , beta1( 1 ) , LipCnst( 0 ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~PrimalDual() override = default;  ///< destructor: it is empty

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 /// (re)initialize the rule in both its roles
 void format( void ) override;

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR DR AND SR COMPUTATION --------------------*/
/*--------------------------------------------------------------------------*/

 void new_deflection( void ) override;

 void new_step( void ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 [[nodiscard]] FunctionValue get_deflection( void ) const override {
  return( vk / tmp_delta );
  }

 /// PrimalDual is its own Deflection
 [[nodiscard]] Deflection * deflection_role( void ) override {
  return( this );
  }

 /// PrimalDual also chooses the stepsize
 [[nodiscard]] bool is_also_stepsize( void ) const override {
  return( true );
  }

 /// the stability center never moves: the steps are from the first point
 [[nodiscard]] bool do_SS( void ) override { return( false ); }

 /// returns the stepsize; the incremental variant is not supported
 FunctionValue get_stepsize( bool incremental = false ) override;

 /// a PrimalDual has no target level: throws
 [[nodiscard]] FunctionValue get_level( void ) const override;

 /// a PrimalDual has no factor beta: throws
 [[nodiscard]] FunctionValue get_beta( void ) const override;

 /// a PrimalDual has no factor beta: throws
 void set_max_beta( FunctionValue alpha ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void update_gamma( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 FunctionValue delta_;     ///< \f$ \Delta_{i-1} \f$
 FunctionValue tmp_delta;  ///< \f$ \Delta_i \f$
 FunctionValue vk;         ///< \f$ \upsilon_i \f$
 FunctionValue step_;      ///< the stepsize \f$ \nu_i \f$
 FunctionValue gamma;      ///< \f$ \gamma \f$
 FunctionValue beta1;      ///< \f$ \hat{\beta}_i \f$
 FunctionValue LipCnst;    ///< the Lipschitz constant

/*--------------------------------------------------------------------------*/
/*------------------------------- FACTORY ----------------------------------*/
/*--------------------------------------------------------------------------*/
// the factory macros cannot be used, as PrimalDual is registered in two
// factories: the same is done by hand

 static class _init { public: _init(); } _initializer;

 [[nodiscard]] const std::string & private_name( void ) const override;

 static const std::string & _private_name( void );

/*--------------------------------------------------------------------------*/

 };  // end( class PrimalDual )

/** @} end( group( PrimalDual_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* PrimalDual.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File PrimalDual.h ---------------------------*/
/*--------------------------------------------------------------------------*/
