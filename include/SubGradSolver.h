/*--------------------------------------------------------------------------*/
/*------------------------ File SubGradSolver.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SubGradSolver class, which implements the CDASolver
 * interface [see CDASolver.h] for the minimization of a convex (or the
 * maximization of a concave) nondifferentiable function, given as a sum of
 * C05Function, by means of a (deflected, projected, incremental)
 * subgradient method, as described in
 *
 *  A. Frangioni, E. Gorgone, B. Gendron. "On the computational efficiency
 *  of subgradient methods: a case study with Lagrangian bounds"
 *  Mathematical Programming Computation 9(4), 573-604, 2017
 *
 * The method is parametric over the stepsize and the deflection rules, which
 * are objects of the classes Stepsize [see Stepsize.h] and Deflection [see
 * Deflection.h] selected by name.
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

#ifndef __SubGradSolver
 #define __SubGradSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CDASolver.h"

#include "C05Function.h"

#include "ColVariable.h"

#include "Deflection.h"

#include "Stepsize.h"

#include <chrono>

#include <random>

#include <unordered_map>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SubGradSolver_CLASSES Classes in SubGradSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS SubGradSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/// a subgradient method for the NDO problems of SMS++
/** The SubGradSolver class implements the CDASolver interface for the
 * minimization of a convex function, or the maximization of a concave one,
 * \f[
 *   \min \left\{ \, f( \lambda ) = f_0( \lambda ) + \sum_{k=1}^K
 *                f_k( \lambda ) \,:\, \lambda \in \Lambda \, \right\}
 * \f]
 * where \f$ \Lambda \f$ is a box (possibly unbounded on either side), each
 * \f$ f_k \f$ is a C05Function and \f$ f_0 \f$ is a LinearFunction or a
 * DQuadFunction. The Block is read as BundleSolver does [see
 * BundleSolver.h]: either it has no sub-Block and its Objective is a
 * FRealObjective whose Function is a C05Function, or each of its sub-Block
 * has a FRealObjective whose Function is a C05Function (a component), the
 * Objective of the Block itself being the 0-th component, if any. The
 * ColVariable of the Block can only be constrained by OneVarConstraint
 * (BoxConstraint, LB0Constraint, NNConstraint) and by their own sign and
 * fixing, which define \f$ \Lambda \f$; the components need not all depend
 * on the same ColVariable. Concave functions to be maximized are internally
 * turned into convex ones to be minimized, and all the values and
 * linearizations described below refer to the minimized function.
 *
 * The algorithm is the unified subgradient method (SM) of
 *
 *  A. Frangioni, E. Gorgone, B. Gendron. "On the computational efficiency
 *  of subgradient methods: a case study with Lagrangian bounds"
 *  Mathematical Programming Computation 9(4), 573-604, 2017
 *
 * which employs the recurrence
 * \f[
 *   \breve{\lambda}_{i+1} = \bar{\lambda}_i - \nu_i d_i \quad , \quad
 *   \lambda_{i+1} = {\rm P}_{\Lambda}( \breve{\lambda}_{i+1} )
 * \f]
 * where \f$ {\rm P}_{\Lambda} \f$ is the orthogonal projection on
 * \f$ \Lambda \f$, \f$ \bar{\lambda}_i \f$ is the *stability center*, which
 * only moves to the new point if the Deflection says so (a *serious step*),
 * the stepsize \f$ \nu_i \f$ is given by a Stepsize object and the direction
 * \f[
 *   d_i = \alpha_i g_i + ( 1 - \alpha_i ) d_{i-1}
 * \f]
 * is a convex combination of the latest subgradient \f$ g_i \f$ and of the
 * previous direction, the deflection coefficient \f$ \alpha_i \f$ being
 * given by a Deflection object (\f$ \alpha_i = 1 \f$ if there is none). The
 * rules are selected by name with the parameters strStepsize and
 * strDeflection, among those registered in the respective factories. The
 * subgradient \f$ g_i \f$ is the sum of one linearization of each component
 * (that of the 0-th one being its gradient). The method also has an
 * *incremental* variant, where a sequence of inner iterations each moves
 * along the subgradient of a single component [see dblIncrFactor].
 *
 * Besides the subgradient, the method keeps track of the linearization
 * errors at the stability center of \f$ g_i \f$ and of \f$ d_i \f$, which
 * enter the stopping condition
 * \f[
 *   t^* \| d_i \| + \epsilon_i \leq \varepsilon \max\{ 1 , | f^{rec}_i | \}
 * \f]
 * where \f$ \epsilon_i \f$ is the linearization error of \f$ d_i \f$,
 * \f$ \| d_i \| \f$ is the norm of the direction projected on the tangent
 * cone of \f$ \Lambda \f$ at the stability center, whatever intProjection
 * is, \f$ f^{rec}_i \f$ the best value found so far,
 * \f$ t^* \f$ is dbltStar and \f$ \varepsilon \f$ the accuracy given by
 * dblRelAcc and dblAbsAcc;
 * as in BundleSolver, this declares optimality when no step of length
 * \f$ t^* \f$ can decrease the function by more than the accuracy, but it
 * does not give a lower bound. The upper bound is the best value found,
 * which is valid whenever compute() stops; a lower bound is only known
 * when it is given by the Block [see Block::get_valid_lower_bound()], or
 * when the direction, which is a convex combination of linearizations
 * \f$ a^\top \lambda + c \f$ of the function (unprojected), has a finite
 * minimum on the box,
 * \f[
 *   \underline{f} = \min \{ \, a^\top \lambda + c \,:\, \lambda \in
 *                   \Lambda \, \} = f( \bar{\lambda} ) - \epsilon +
 *   \min \{ \, a^\top ( \lambda - \bar{\lambda} ) \,:\, \lambda \in
 *                \Lambda \, \}
 * \f]
 * which is the bound of BundleSolver on the box of the variables: it is
 * always finite if all the variables have finite bounds, and only by chance
 * otherwise (for a subgradient method, hardly ever), hence normally the lower
 * bound stays - INF. The bound is also \f$ f( \bar{\lambda} ) - \epsilon \f$
 * when the norm of the projected direction is "0" [see dblNZEps]. When
 * either bound is close enough to the best value, optimality is
 * *certified*. When the method stalls, i.e., the stepsize stays below
 * 1e-8 * dbltStar for 100 consecutive iterations, compute() returns
 * kLowPrecision.
 *
 * For each component that has a global pool of linearizations
 * [see C05Function::intGPMaxSz], SubGradSolver keeps there the latest
 * linearization and a convex combination of the past ones, and
 * get_dual_solution() sets the latter as the important linearization of
 * the component [see C05Function::set_important_linearization()]. The
 * combination [see intPrimalRec] is either that of the direction, whose
 * weights are the deflection coefficients, or the ergodic average, whose
 * weights are the stepsizes, which is the only one without a Deflection.
 * When the components are LagBFunction, this is a convex combination of
 * the solutions of the Lagrangian subproblems, which both the Volume-type
 * deflections and the ergodic average are known to make converge to an
 * optimal solution of the convexified relaxation; this is how
 * LagrangianDualSolver recovers its primal solution.
 *
 * Only box constraints on the variables are handled, by projection. */

class SubGradSolver : public CDASolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- FRIEND CLASSES ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Friend classes
 *
 * Stepsize and Deflection are friends of SubGradSolver, so that they can
 * read its data, and they give the classes derived from them access to it
 * through protected methods [see Stepsize.h and Deflection.h].
 *  @{ */

 friend class Stepsize;
 friend class Deflection;

/** @} ---------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *  @{ */

 using Index = Function::Index;                       ///< an index
 using c_Index = Function::c_Index;                   ///< a const index
 using VarValue = Function::FunctionValue;            ///< a value
 using c_VarValue = Function::c_FunctionValue;        ///< a const value
 using Vec_VarValue = Function::Vec_FunctionValue;    ///< a vector of values
 using LinearCombination = C05Function::LinearCombination;
                                         ///< a combination of linearizations

/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters of SubGradSolver

 enum int_par_type_SGS {
  intProjection = intLastParCDAS ,  ///< which vectors are projected
  intScheme ,          ///< safe rule and stepsize / deflection restriction
  intKeepHat ,         ///< whether \f$ \hat{\lambda}_i \f$ is kept
  intSeed ,            ///< seed of the shuffle of the incremental variant
  intResetAlg ,        ///< what is reset at each compute()
  intPrimalRec ,       ///< the combination kept for the primal recovery
  intColorGreen ,      ///< ColorTV: consecutive green iterations
  intColorYellow ,     ///< ColorTV: consecutive yellow iterations
  intColorRed ,        ///< ColorTV: consecutive red iterations
  intFumeroR1 ,        ///< FumeroTV: the parameter \f$ r_1 \f$
  intFumeroEta1 ,      ///< FumeroTV: failures in phase two
  intFumeroEta2 ,      ///< FumeroTV: failures in phase one
  intVolumeTauIter ,   ///< Volume: iterations between checks of tau
  intPDAverage ,       ///< PrimalDual: simple (1) or weighted (0) averages
  intLastParSGS        ///< first allowed parameter value for derived classes
  };

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters of SubGradSolver

 enum dbl_par_type_SGS {
  dbltStar = dblLastParCDAS ,  ///< the longest step "that matters"
  dblNZEps ,            ///< when the projected direction is "0"
  dblIncrFactor ,       ///< the incremental factor
  dblLipFactor ,        ///< the scaling factor of the Lipschitz constant
  dblPolyakBeta ,       ///< Polyak: the factor \f$ \beta \f$
  dblColorBeta0 ,       ///< ColorTV: the initial \f$ \beta \f$
  dblFumeroBeta0 ,      ///< FumeroTV: the initial \f$ \beta \f$
  dblFumeroSigmaMin ,   ///< FumeroTV: the threshold of \f$ \sigma( r ) \f$
  dblVolumeTauInit ,    ///< Volume: the initial \f$ \tau \f$
  dblVolumeTauMin ,     ///< Volume: the minimum \f$ \tau \f$
  dblVolumeTauFactor ,  ///< Volume: the decrease factor of \f$ \tau \f$
  dblVolumeM ,          ///< Volume: the descent parameter
  dblLastParSGS         ///< first allowed parameter value for derived classes
  };

/*--------------------------------------------------------------------------*/
 /// public enum for the string algorithmic parameters of SubGradSolver

 enum str_par_type_SGS {
  strStepsize = strLastParCDAS ,  ///< name of the Stepsize
  strDeflection ,                 ///< name of the Deflection
  strLastParSGS   ///< first allowed parameter value for derived classes
  };

/** @} ---------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing SubGradSolver
 *  @{ */

 /// constructor: all the parameters are given their default value
 SubGradSolver( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 /// destructor: cleanly detaches the SubGradSolver from the Block
 ~SubGradSolver() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set the Block that the SubGradSolver has to solve
 /** Set the Block that the SubGradSolver has to solve, which must have the
  * structure described in the general notes of the class; throws
  * std::invalid_argument or std::logic_error if it has not. The initial
  * stability center is the current value of the ColVariable, projected on
  * \f$ \Lambda \f$. The global pool of every component that can have one
  * [see C05Function::intGPMaxSz] is enlarged by 3 linearizations, which are
  * those SubGradSolver uses, so that it does not touch those of any other
  * Solver using the same component. */

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 using ThinComputeInterface::set_par;  // restore the hidden overloads
 using ThinComputeInterface::get_int_par;
 using ThinComputeInterface::get_dbl_par;
 using ThinComputeInterface::get_str_par;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the int parameters of SubGradSolver
 /** Set the int parameters specific of SubGradSolver, together with those
  * of CDASolver that it "listens to":
  *
  * - intMaxIter [Inf< int >()]: the maximum number of iterations of each
  *   compute(), inner iterations of the incremental variant included
  *
  * - intLogVerb [0]: the verbosity of the log, 0 = none, 1 = the final
  *   state of each compute(), 2 = one line per iteration, 3 = details of
  *   each iteration
  *
  * - intEverykIt [0]: how often the eEverykIteration events are called
  *
  * - intProjection [0]: which of the vectors \f$ g_i \f$, \f$ d_{i-1} \f$
  *   and \f$ d_i \f$ are projected on the tangent cone of \f$ \Lambda \f$ at
  *   the stability center, coded bit-wise: bit 0 for \f$ g_i \f$, bit 1 for
  *   \f$ d_{i-1} \f$, bit 2 for \f$ d_i \f$; asking both \f$ g_i \f$ and
  *   \f$ d_{i-1} \f$ makes asking \f$ d_i \f$ useless, and without a
  *   Deflection only bit 0 is considered
  *
  * - intScheme [1]: the convergence scheme, coded bit-wise: bit 0 says
  *   whether the *safe rule* is used, bit 1 whether the scheme is
  *   *stepsize-restricted* (the stepsize is computed after the direction,
  *   and the safe rule bounds it by the deflection coefficient) or
  *   *deflection-restricted* (the stepsize is computed first, and the safe
  *   rule bounds the deflection coefficient by it); without a Deflection
  *   the scheme is always stepsize-restricted
  *
  * - intKeepHat [1]: whether the point \f$ \hat{\lambda}_{i+1} =
  *   \alpha_{i+1} \lambda_i + ( 1 - \alpha_{i+1} ) \hat{\lambda}_i \f$ is
  *   kept, whose linearization error is then also used in the stopping
  *   condition; it is ignored by the incremental variant
  *
  * - intSeed [0]: the seed of the random shuffle of the components in the
  *   incremental variant
  *
  * - intResetAlg [0]: what is reset at the beginning of each compute()
  *   after the first, coded bit-wise: bit 0 for the rules (as after a
  *   change of the function), bit 1 for the stability center, which is
  *   then read from the value of the ColVariable
  *
  * - intPrimalRec [0]: the convex combination of linearizations that is
  *   kept in the global pools for the primal recovery [see
  *   get_dual_solution()]: 0 for that of the direction, whose weights are
  *   the deflection coefficients (as in the Volume algorithm and in the
  *   Primal-Dual one), 1 for the ergodic average, whose weights are the
  *   stepsizes; without a Deflection the ergodic average is the only one,
  *   whatever the value
  *
  * - intColorGreen [1], intColorYellow [400], intColorRed [10]: the
  *   number of consecutive iterations of each color after which ColorTV
  *   changes \f$ \beta \f$ [see ColorTV.h]
  *
  * - intFumeroR1 [5], intFumeroEta1 [5], intFumeroEta2 [5]: the parameters
  *   \f$ r_1 \f$, \f$ \eta_1 \f$ and \f$ \eta_2 \f$ of FumeroTV [see
  *   FumeroTV.h]
  *
  * - intVolumeTauIter [100]: the iterations between two checks of
  *   \f$ \tau \f$ in Volume [see Volume.h]
  *
  * - intPDAverage [1]: simple (1) or weighted (0) averages in PrimalDual
  *   [see PrimalDual.h] */

 void set_par( idx_type par , int value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the double parameters of SubGradSolver
 /** Set the double parameters specific of SubGradSolver, together with
  * those of CDASolver that it "listens to":
  *
  * - dblMaxTime [Inf< double >()]: the maximum time of each compute(), in
  *   seconds
  *
  * - dblEveryTTm [0]: how often the eEveryTTime events are called
  *
  * - dblRelAcc [1e-6]: the relative accuracy of the stopping conditions
  *
  * - dblAbsAcc [Inf< double >()]: the absolute accuracy of the stopping
  *   conditions, the smallest of the two being used
  *
  * - dblUpCutOff [Inf< double >()], dblLwCutOff [- Inf< double >()]: the
  *   cutoffs [see Solver.h], which stop compute() with kCutOff
  *
  * - dbltStar [100]: the longest step "that matters" in the stopping
  *   condition [see the general notes of the class]
  *
  * - dblNZEps [0]: when the norm of the direction, projected on the tangent
  *   cone of the box at the stability center, is at most dblNZEps, it is
  *   taken as 0, and the value at the center minus the linearization error
  *   is declared a lower bound; as in BundleSolver, the default 0 only does
  *   it when the projection is exactly 0, while a too large value may give
  *   a lower bound which is not valid
  *
  * - dblIncrFactor [0]: the incremental factor; after each full (outer)
  *   iteration, ceil( dblIncrFactor * K ) inner iterations are done along
  *   the subgradient of one component only, K being the number of
  *   components (the 0-th one included); the incremental variant cannot be
  *   deflected
  *
  * - dblLipFactor [1]: the scaling factor \f$ \chi \f$ of the Lipschitz
  *   constant in the stepsize of the incremental variant, and the factor F
  *   of PrimalDual
  *
  * - dblPolyakBeta [1e-2]: the factor \f$ \beta \f$ of Polyak
  *
  * - dblColorBeta0 [1]: the initial \f$ \beta \f$ of ColorTV
  *
  * - dblFumeroBeta0 [1], dblFumeroSigmaMin [1e-4]: the initial \f$ \beta \f$
  *   and the threshold of \f$ \sigma( r ) \f$ of FumeroTV
  *
  * - dblVolumeTauInit [1], dblVolumeTauMin [1e-4], dblVolumeTauFactor
  *   [0.5], dblVolumeM [0.1]: the parameters of Volume */

 void set_par( idx_type par , double value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the string parameters of SubGradSolver
 /** Set the string parameters specific of SubGradSolver:
  *
  * - strStepsize ["ColorTV"]: the name of the Stepsize class, among Polyak,
  *   ColorTV, FumeroTV and PrimalDual (or any other registered in the
  *   Stepsize factory)
  *
  * - strDeflection ["Volume"]: the name of the Deflection class, among
  *   Volume and PrimalDual (or any other registered in the Deflection
  *   factory); "" or "none" for no deflection. A rule that is both a
  *   Stepsize and a Deflection, as PrimalDual, plays both roles when named
  *   in strStepsize, and then strDeflection must either be empty or name
  *   the same class
  *
  * The rules are constructed by set_Block(), or by the next compute() if
  * they are changed afterwards. */

 void set_par( idx_type par , std::string && value ) override;

/*--------------------------------------------------------------------------*/
 /// set the ostream for the SubGradSolver log
 void set_log( std::ostream * log_stream = nullptr ) override {
  f_log = log_stream;
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem encoded by the current Block
 *  @{ */

 /// (try to) solve the problem encoded in the Block
 /** Runs the subgradient method from the current stability center. The
  * return value is
  *
  * - kOK if the stopping condition, or one of the certified ones, is met
  *   [see the general notes of the class];
  *
  * - kUnbounded if the function is unbounded below (above for a concave
  *   one), which is detected either by a component evaluating at - INF or
  *   by the best value going below the conditional valid lower bound of the
  *   Block [see Block::get_valid_lower_bound()];
  *
  * - kStopIter / kStopTime if intMaxIter / dblMaxTime is reached;
  *
  * - kCutOff if one of the cutoffs is reached;
  *
  * - kLowPrecision if the method stalls, i.e., the stepsize stays tiny;
  *
  * - kError if a component cannot be evaluated, or it only gives vertical
  *   linearizations, which SubGradSolver cannot handle;
  *
  * - kBlockLocked if the Block cannot be read_lock()-ed to process the
  *   outstanding Modification. */

 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 *  @{ */

 /// returns a valid lower bound on the optimal value
 /** For a convex function to be minimized, this is the best valid lower
  * bound known [see the general notes of the class], - INF if none, which
  * is the usual case for a subgradient method; for a concave one to be
  * maximized, it is the best value found, i.e., that of the point that
  * get_var_solution() writes, which is valid whenever compute() stops. */

 OFValue get_lb( void ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns a valid upper bound on the optimal value
 /** The symmetric of get_lb(). */

 OFValue get_ub( void ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if a point has been evaluated

 bool has_var_solution( void ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the value of the point has been finitely evaluated

 bool is_var_feasible( void ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the point found in the ColVariable
 /** Writes in the ColVariable the best point found (the stability center
  * if no point has been evaluated yet). \p solc is ignored. */

 void get_var_solution( Configuration * solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the dual solution can be given
 /** The dual solution is, for each component, the convex combination of
  * linearizations in its global pool chosen by intPrimalRec (or the latest
  * linearization if the combination has been lost): it can be given if
  * every component has a global pool and that linearization is still
  * there. */

 bool has_dual_solution( void ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the certified lower bound is finite

 bool is_dual_feasible( void ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the dual solution in the components
 /** Sets, for each component, the convex combination of linearizations in
  * its global pool chosen by intPrimalRec (or the latest linearization if
  * the combination has been lost) as its important linearization [see
  * C05Function::set_important_linearization()]. If \p solc is a
  * SimpleConfiguration< std::vector< int > >, only the components whose
  * indices are listed there are considered. Throws std::logic_error if a
  * component has no global pool. */

 void get_dual_solution( Configuration * solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the elapsed time since the start of the last compute()

 double get_elapsed_time( void ) const override {
  std::chrono::duration< double > elapsed =
   std::chrono::system_clock::now() - c_start;
  return( elapsed.count() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the iterations of the last compute()

 long get_elapsed_iterations( void ) const override { return( ParIter ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the number of calls to compute()

 long get_elapsed_calls( void ) const override { return( SCalls ); }

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR EVENTS HANDLING -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Set event handlers
 *
 * SubGradSolver supports the events eBeforeTermination (called just before
 * an optimality stop), eEverykIteration and eEveryTTime.
 *  @{ */

 [[nodiscard]] EventID max_event_number( void ) const override {
  return( 3 );
  }

/** @} ---------------------------------------------------------------------*/
/*----------------- METHODS FOR ACCESSING THE DATA OF THE Block ------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the data of the Block
 *  @{ */

 /// returns the number of components, the 0-th one excluded
 [[nodiscard]] Index n_components( void ) const { return( NrFi ); }

 /// returns the i-th component
 [[nodiscard]] C05Function * component( Index i ) const {
  return( v_c05f[ i ] );
  }

 /// returns the 0-th component (nullptr if there is none)
 [[nodiscard]] C05Function * zeroth_component( void ) const {
  return( f_zero );
  }

 /// returns the number of ColVariable the components depend upon
 [[nodiscard]] Index get_num_var( void ) const { return( NumVar ); }

 /// returns the current Stepsize, if any
 [[nodiscard]] Stepsize * get_Stepsize( void ) const { return( stepsize ); }

 /// returns the current Deflection, if any
 [[nodiscard]] Deflection * get_Deflection( void ) const {
  return( deflection );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the value at the stability center (+ INF if not known)
 [[nodiscard]] VarValue get_Fi_bar( void ) const { return( FiBar ); }

 /// the best value found so far (+ INF if none)
 [[nodiscard]] VarValue get_Fi_best( void ) const { return( FiBest ); }

 /// the value at \f$ \hat{\lambda}_i \f$ (+ INF if not kept)
 [[nodiscard]] VarValue get_Fi_hat( void ) const {
  return( SGPar4 ? FiHat : Inf< VarValue >() );
  }

 /// the iterations since the rules were last formatted
 [[nodiscard]] Index get_tot_iter( void ) const { return( TotIter ); }

/** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of SubGradSolver
 *  @{ */

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastParSGS ) );
  }

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( idx_type( dblLastParSGS ) );
  }

 [[nodiscard]] idx_type get_num_str_par( void ) const override {
  return( idx_type( strLastParSGS ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 [[nodiscard]] const std::string & get_str_par( idx_type par )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override;

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods
 *  @{ */

 /// reads the structure of the Block: components, ColVariable, bounds
 void build_structure( bool generate = true );

 /// forgets the structure of the Block
 void clear_structure( void );

 /// (re)reads the bounds of the ColVariable
 void read_bounds( void );

 /// processes the outstanding Modification, returns true if all changes
 bool process_outstanding_Modification( void );

 /// (re)constructs the rules out of strStepsize and strDeflection
 void build_rules( void );

 /// deletes the rules
 void delete_rules( void );

 /// forgets all that is known about the function, keeps the center
 void reset_function_info( void );

 /// re-initializes the rules and the direction
 void init_deflection( void );

/*--------------------------------------------------------------------------*/

 /// evaluates the function (or component wFi) at Lambda
 bool Fi_and_Gi( Index wFi = Inf< Index >() );

 /// evaluates the h-th function (0-th component last), adds to Gi
 bool eval_function( Index h , VarValue & val , VarValue & cnst );

 /// computes the deflection coefficient [and the stepsize]
 void form_D( void );

 /// computes the direction [and the stepsize]
 void save_dir( void );

 /// moves the stability center to Lambda
 void goto_lambda1( void );

 /// computes the (unprojected) tentative point
 void form_lambda1( void );

 /// projects a point on the box, returns true if it has changed
 bool proj_fsb( Vec_VarValue & Lpt );

 /// projects a vector on the tangent cone of the box at Lpt
 void proj_tng( Vec_VarValue & Gpt , const Vec_VarValue & Lpt );

 /// updates Sigma at a null step
 void update_sigma( void );

 /// updates LambdaHat
 void chg_lambda_hat( void );

 /// updates the lower bounds known by the Block and the components
 void updt_lower_bound( void );

 /// updates the certified lower bound out of the aggregate linearization
 void update_cert_bound( void );

 /// the heuristic stopping condition
 [[nodiscard]] bool is_optimal( void ) const;

 /// the accuracy for the value v
 [[nodiscard]] VarValue max_error( VarValue v ) const {
  return( std::min( RelAcc * std::max( VarValue( 1 ) , std::abs( v ) ) ,
		    AbsAcc ) );
  }

 /// the best true lower bound on the optimal value
 [[nodiscard]] VarValue global_lb( void ) const {
  return( std::max( TrueLB , CertLB ) );
  }

 /// the value of the point given by get_var_solution()
 [[nodiscard]] VarValue sol_value( void ) const {
  return( FiBest );
  }

 /// the (estimate of the) global Lipschitz constant
 [[nodiscard]] VarValue lipschitz( void ) const;

 /// stores the combination for the primal recovery in the global pools
 void store_primal( void );

 /// writes Lambda in the ColVariable
 void lambda_to_vars( const Vec_VarValue & L );

 /// convex-minimization sign of a value of the original function
 [[nodiscard]] VarValue rs( VarValue v ) const {
  return( f_convex ? v : - v );
  }

/** @} ---------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - -

 int MaxIter;          ///< maximum number of iterations
 int LogVerb;          ///< verbosity of the log
 int EverykIt;         ///< periodicity of eEverykIteration events
 int SGPar1;           ///< the vectors to be projected
 int SGPar3;           ///< the convergence scheme
 int SGPar4;           ///< whether LambdaHat is kept
 int SGPar5;           ///< the seed
 int RstAlg;           ///< what is reset at each compute()
 int PrimalRec;        ///< the combination kept for the primal recovery

 double MaxTime;       ///< maximum time of each compute()
 double EveryTTm;      ///< periodicity of eEveryTTime events
 double RelAcc;        ///< relative accuracy
 double AbsAcc;        ///< absolute accuracy
 double UpCutOff;      ///< upper cutoff
 double LwCutOff;      ///< lower cutoff
 double tStar;         ///< the longest step that matters
 double NZEps;         ///< when the projected direction is "0"
 double SGPar2;        ///< the incremental factor

 std::string StepsizeName;    ///< name of the Stepsize
 std::string DeflectionName;  ///< name of the Deflection

 /// the int parameters of the rules, from intColorGreen on
 std::vector< int > v_rule_int;

 /// the double parameters of the rules, from dblLipFactor on
 std::vector< double > v_rule_dbl;

 // the structure of the Block- - - - - - - - - - - - - - - - - - - - - - - -

 bool f_convex;                       ///< true if the function is convex
 std::vector< C05Function * > v_c05f;  ///< the components
 C05Function * f_zero;                ///< the 0-th component, if any
 Index NrFi;                          ///< number of components

 std::vector< ColVariable * > LamVcblr;  ///< the ColVariable
 Index NumVar;                           ///< their number

 /// for each function (components, then the 0-th) the position of each of
 /// its active Variable in LamVcblr; empty if it is the identity
 std::vector< std::vector< Index > > v_loc2glob;

 Vec_VarValue lb;  ///< the lower bounds on the ColVariable
 Vec_VarValue ub;  ///< the upper bounds on the ColVariable

 std::vector< bool > v_pool;    ///< whether component k has a global pool
 std::vector< bool > v_dir_ok;  ///< whether the pool of k holds the direction
 std::vector< Index > v_slot;   ///< the first name in the pool of k

 /// the first name in the pool of each C05Function SubGradSolver has used,
 /// kept across changes of the structure (the pointers are only keys)
 std::unordered_map< const C05Function * , Index > f_slot_base;

 bool f_empty_box;  ///< true if some lower bound is larger than its upper

 Vec_VarValue v_buf;  ///< buffer for the linearizations of one function

 // the rules- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Stepsize * stepsize;      ///< the Stepsize
 Deflection * deflection;  ///< the Deflection, nullptr if none
 bool f_shared_rule;       ///< true if they are the same object
 bool f_rules_dirty;       ///< true if the rules have to be reconstructed
 int f_scheme;             ///< the convergence scheme actually used

 // the state of the algorithm- - - - - - - - - - - - - - - - - - - - - - - -

 int Result;           ///< result of the latest call to compute()
 Index SCalls;         ///< number of calls to compute()
 Index ParIter;        ///< iterations in the current call to compute()
 Index TotIter;        ///< iterations since the rules were formatted
 Index FiEvaltns;      ///< total number of evaluations of the function
 Index GiEvaltns;      ///< total number of subgradients

 Vec_VarValue LambdaBar;   ///< the stability center
 Vec_VarValue Lambda;      ///< the current (tentative) point
 Vec_VarValue LambdaBest;  ///< the best point found so far
 Vec_VarValue LambdaHat;   ///< the point \f$ \hat{\lambda}_i \f$

 VarValue FiBar;       ///< the value at LambdaBar
 VarValue FiLambda;    ///< the value at Lambda
 VarValue FiBest;      ///< the best value found so far
 VarValue FiHat;       ///< the value at LambdaHat

 Vec_VarValue Gi;      ///< the subgradient (possibly projected)
 Vec_VarValue GiRaw;   ///< the subgradient, never projected
 VarValue GiCnst;      ///< the constant of the linearization GiRaw
 Vec_VarValue dir;     ///< the direction
 Vec_VarValue dirRaw;  ///< the aggregate of the GiRaw
 VarValue dirCnst;     ///< the constant of the linearization dirRaw

 VarValue alpha;       ///< the deflection coefficient
 VarValue step;        ///< the stepsize

 VarValue Sigma;       ///< linearization error of Gi at LambdaBar
 VarValue Epsilon;     ///< linearization error of dir at LambdaBar
 VarValue SigmaHat;    ///< linearization error of Gi at LambdaHat
 VarValue HatEpsilon;  ///< linearization error of dir at LambdaHat

 VarValue NrmGi;       ///< the squared norm of Gi
 VarValue NrmDir;      ///< the squared norm of dir
 VarValue NrmPrj;      ///< the squared norm of the projected direction
 VarValue dGk;         ///< the scalar product of dir and Gi
 VarValue dM1Gk;       ///< the scalar product of the previous dir and Gi
 bool dM1GkDone;       ///< true if dM1Gk has been computed
 VarValue MaxNrmGi;    ///< the largest norm of a full subgradient

 bool LHasProj;        ///< true if Lambda has been projected
 bool DoSS;            ///< true if the last step was a serious one
 bool f_dir_valid;     ///< true if dir is a valid aggregate
 VarValue ErgSum;      ///< the sum of the stepsizes of the ergodic average
 bool f_eval_ok;       ///< true if a point has been evaluated

 Index DirPos;         ///< position of the direction in the global pools

 VarValue TrueLB;      ///< the best true lower bound from the Block
 VarValue CondLB;      ///< the best conditional lower bound
 VarValue CertLB;      ///< the certified lower bound

 Index CSSCntr;        ///< consecutive serious steps
 Index CNSCntr;        ///< consecutive null steps
 Index CSmallStep;     ///< consecutive small steps

 Index NItIncr;        ///< inner iterations after an outer one
 bool InnIter;         ///< true if the current iteration is an inner one
 std::vector< Index > Seq;    ///< the shuffled functions
 std::mt19937 myrandom;       ///< the random generator for the shuffle

 bool f_need_rebuild;  ///< true if the structure has to be read again

 std::chrono::time_point< std::chrono::system_clock > c_start;
                                        ///< the start of the last compute()

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class SubGradSolver )

/** @} end( group( SubGradSolver_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SubGradSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File SubGradSolver.h --------------------------*/
/*--------------------------------------------------------------------------*/
