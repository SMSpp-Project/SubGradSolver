/*--------------------------------------------------------------------------*/
/*--------------------------- File Volume.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the class Volume, which implements the Deflection interface
 * [see Deflection.h] with the deflection rule of the modified Volume
 * algorithm, i.e., that of the "poorman's" bundle method. For details see
 *
 *  L. Bahiense, N. Maculan, C. Sagastizabal. "The volume algorithm
 *  revisited: relation with bundle methods"
 *  Mathematical Programming 94, 41-69, 2002
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

#ifndef __Volume
 #define __Volume  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Deflection.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Volume_CLASSES Classes in Volume.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS Volume --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the deflection rule of the modified Volume algorithm
/** The class Volume implements the deflection rule of the modified Volume
 * algorithm. The deflection coefficient solves the one-dimensional
 * quadratic problem
 * \f[
 *   \alpha_i = \arg\min \left\{ \, \nu_{i-1} \| \alpha g_i
 *                    + ( 1 - \alpha ) d_{i-1} \|^2 / 2
 *                    + \alpha \sigma_i + ( 1 - \alpha ) \epsilon_{i-1}
 *                    \,:\, \alpha \in [ 0 , 1 ] \, \right\}
 * \f]
 * where \f$ \sigma_i \f$ and \f$ \epsilon_{i-1} \f$ are the linearization
 * errors at the stability center of \f$ g_i \f$ and \f$ d_{i-1} \f$, and
 * \f$ \nu_{i-1} \f$ the last stepsize; the solution is safeguarded in
 * \f$ [ \tau / 10 , \tau ] \f$, where \f$ \tau \f$ starts at dblVolumeTauInit
 * and is multiplied by dblVolumeTauFactor (not below dblVolumeTauMin) every
 * intVolumeTauIter iterations in which the best value has not decreased by
 * at least 1%. A serious step takes place when the function decreases by at
 * least dblVolumeM times dblRelAcc times the best value [see delta()]. */

class Volume : public Deflection
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------- CONSTRUCTOR AND DESTRUCTOR ------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor: takes the SubGradSolver using the rule
 explicit Volume( SubGradSolver * slvr = nullptr )
  : Deflection( slvr ) , tau( 1 ) , lastvalue( Inf< FunctionValue >() ) ,
    alpha( 1 ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~Volume() override = default;  ///< destructor: it is empty

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 void format( void ) override;

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR DEFLECTION COMPUTATION ------------------*/
/*--------------------------------------------------------------------------*/

 void new_deflection( void ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 [[nodiscard]] FunctionValue get_deflection( void ) const override {
  return( alpha );
  }

 [[nodiscard]] bool do_SS( void ) override;

 [[nodiscard]] FunctionValue delta( void ) const override;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 FunctionValue tau;        ///< the safeguard \f$ \tau \f$
 FunctionValue lastvalue;  ///< the best value when tau was last checked
 FunctionValue alpha;      ///< the deflection coefficient

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class Volume )

/** @} end( group( Volume_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Volume.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File Volume.h ------------------------------*/
/*--------------------------------------------------------------------------*/
