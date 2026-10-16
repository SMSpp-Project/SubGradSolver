/*--------------------------------------------------------------------------*/
/*--------------------------- File Polyak.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the class Polyak, which implements the Stepsize interface
 * [see Stepsize.h] with the simplest target value stepsize rule, whose
 * factor \f$ \beta \f$ is constant and whose level is the best known lower
 * bound on the function. For details see
 *
 *  B.T. Polyak. "Minimization of unsmooth functionals"
 *  Zh. Vychisl. Mat. Fiz. 9(3), 509-521, 1969
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

#ifndef __Polyak
 #define __Polyak  /* self-identification: #endif at the end of the file */

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
/** @defgroup Polyak_CLASSES Classes in Polyak.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS Polyak --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Polyak target value stepsize rule
/** The class Polyak implements a target value stepsize rule whose
 * \f$ \beta_i = \beta \f$ is constant (the parameter dblPolyakBeta of
 * SubGradSolver) and whose level \f$ f^{lev}_i \f$ is the best known lower
 * bound on the function, which therefore must be finite: the bound is the
 * (possibly conditional) valid bound of the Block, or the sum of those of
 * the components [see SubGradSolver]. The rule is the one with the
 * strongest convergence properties when the bound is the optimal value, and
 * it can be quite slow when it is not tight. */

class Polyak : public Stepsize
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------- CONSTRUCTOR AND DESTRUCTOR ------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor: takes the SubGradSolver using the rule
 explicit Polyak( SubGradSolver * slvr = nullptr ) : Stepsize( slvr ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~Polyak() override = default;  ///< destructor: it is empty

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

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class Polyak )

/** @} end( group( Polyak_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Polyak.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File Polyak.h ------------------------------*/
/*--------------------------------------------------------------------------*/
