/*--------------------------------------------------------------------------*/
/*--------------------------- File ColorTV.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the class ColorTV, which implements the Stepsize interface
 * [see Stepsize.h] with the target value stepsize rule of the original
 * Volume algorithm, which classifies the iterations as "green", "yellow" or
 * "red" according to the improvement of the function and changes
 * \f$ \beta_i \f$ accordingly. For details see
 *
 *  F. Barahona, R. Anbil. "The volume algorithm: producing primal solutions
 *  with a subgradient method"
 *  Mathematical Programming 87(3), 385-399, 2000
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

#ifndef __ColorTV
 #define __ColorTV  /* self-identification: #endif at the end of the file */

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
/** @defgroup ColorTV_CLASSES Classes in ColorTV.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS ColorTV -------------------------------*/
/*--------------------------------------------------------------------------*/
/// the target value stepsize rule of the Volume algorithm
/** The class ColorTV implements the target value stepsize rule of the
 * Volume algorithm. Each iteration is "green" if \f$ d_i^\top g_i > 0 \f$
 * and the function has decreased by a relative amount of at least 1e-6,
 * "yellow" if \f$ d_i^\top g_i \leq 0 \f$ but the function has decreased,
 * "red" otherwise; \f$ \beta_i \f$ (starting from dblColorBeta0) is doubled
 * after intColorGreen consecutive green iterations (as long as it is < 2),
 * multiplied by 1.1 after intColorYellow consecutive yellow ones, and
 * multiplied by 0.67 after intColorRed consecutive red ones (as long as it
 * is > 5e-4). The level \f$ f^{lev}_i \f$ is kept 5% below the current
 * value, but never below the best known lower bound on the function. Unlike
 * the other target value rules, ColorTV also works when no lower bound is
 * known. */

class ColorTV : public Stepsize
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// the colors of an iteration
 enum condition { green ,   ///< good step
		  yellow ,  ///< neither good nor bad step
		  red       ///< bad step
                  };

/*--------------------------------------------------------------------------*/
/*---------------------- CONSTRUCTOR AND DESTRUCTOR ------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor: takes the SubGradSolver using the rule
 explicit ColorTV( SubGradSolver * slvr = nullptr )
  : Stepsize( slvr ) , lastgreeniter( 0 ) , lastyellowiter( 0 ) ,
    lastrediter( 0 ) , NrIter( 0 ) , LwrBnd( - Inf< FunctionValue >() ) ,
    lastvalue( Inf< FunctionValue >() ) ,
    FiLmbd( Inf< FunctionValue >() ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~ColorTV() override = default;  ///< destructor: it is empty

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 void format( void ) override;

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
/*--------------------------------------------------------------------------*/

 void new_step( void ) override;

 [[nodiscard]] bool needs_dM1Gk( void ) const override { return( true ); }

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

 condition coloring( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 Index lastgreeniter;       ///< the last green iteration
 Index lastyellowiter;      ///< the last yellow iteration
 Index lastrediter;         ///< the last red iteration
 Index NrIter;              ///< the iterations since the last format()

 FunctionValue LwrBnd;      ///< the best known lower bound
 FunctionValue lastvalue;   ///< the value at the center at the last step
 FunctionValue FiLmbd;      ///< the value at the last evaluated point

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class ColorTV )

/** @} end( group( ColorTV_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ColorTV.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File ColorTV.h -----------------------------*/
/*--------------------------------------------------------------------------*/
