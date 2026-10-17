/*--------------------------------------------------------------------------*/
/*--------------------------- File FumeroTV.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the class FumeroTV, which implements the Stepsize interface
 * [see Stepsize.h] with a target value stepsize rule whose level starts at
 * the best known lower bound and moves exponentially towards the best value
 * found so far. For details see
 *
 *  F. Fumero. "A modified subgradient algorithm for Lagrangean relaxation"
 *  Computers and Operations Research 28(1), 33-52, 2001
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

#ifndef __FumeroTV
 #define __FumeroTV  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

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
/** @defgroup FumeroTV_CLASSES Classes in FumeroTV.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS FumeroTV -------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Fumero target value stepsize rule
/** The class FumeroTV implements a target value stepsize rule whose level
 * is
 * \f[
 *   f^{lev}_i = \sigma( r ) \underline{f} + ( 1 - \sigma( r ) ) f^{rec}_i
 *   \quad , \quad
 *   \sigma( r ) = e^{ - 0.6933 ( r / r_1 )^{3.26} }
 * \f]
 * where \f$ \underline{f} \f$ is the best known lower bound on the function
 * (which therefore must be finite), \f$ f^{rec}_i \f$ the best value found
 * so far and \f$ r \f$ counts the failures to improve it. The method works
 * in two phases: in the first one, after intFumeroEta2 consecutive failures
 * \f$ \beta \f$ is decreased as \f$ \beta / ( 2 \beta + 1 ) \f$ and \f$ r
 * \f$ increased, until \f$ \sigma( r ) \f$ falls below dblFumeroSigmaMin;
 * in the second one \f$ \beta \f$ is halved after intFumeroEta1 consecutive
 * failures and doubled at each improvement. The initial \f$ \beta \f$ is
 * dblFumeroBeta0, \f$ r_1 \f$ is intFumeroR1; the names of the parameters
 * follow \f$ \bar{\eta}_1 \f$ and \f$ \bar{\eta}_2 \f$ of the paper. */

class FumeroTV : public Stepsize
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------- CONSTRUCTOR AND DESTRUCTOR ------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor: takes the SubGradSolver using the rule
 explicit FumeroTV( SubGradSolver * slvr = nullptr )
  : Stepsize( slvr ) , FiLambda( Inf< FunctionValue >() ) ,
    LwrBnd( - Inf< FunctionValue >() ) , FiRef( Inf< FunctionValue >() ) ,
    sigma( 1 ) , rc( 0 ) , r2( 0 ) , etac( 0 ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~FumeroTV() override = default;  ///< destructor: it is empty

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 void format( void ) override;

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
/*--------------------------------------------------------------------------*/

 /// sets \f$ \beta \f$ and the level; throws if no lower bound is known
 void new_step( void ) override;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 bool update_target_level( void ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 [[nodiscard]] FunctionValue exp_fun( Index ri ) const;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 FunctionValue FiLambda;   ///< the value at the last evaluated point
 FunctionValue LwrBnd;     ///< the best known lower bound
 FunctionValue FiRef;      ///< the best value found so far
 FunctionValue sigma;      ///< the current value of \f$ \sigma( r ) \f$

 Index rc;                 ///< the current value of \f$ r \f$
 Index r2;                 ///< the value of \f$ r \f$ ending phase one
 Index etac;               ///< the consecutive failures

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class FumeroTV )

/** @} end( group( FumeroTV_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* FumeroTV.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File FumeroTV.h -----------------------------*/
/*--------------------------------------------------------------------------*/
