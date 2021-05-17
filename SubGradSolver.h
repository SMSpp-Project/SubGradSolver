/*--------------------------------------------------------------------------*/
/*---------------------- File SubGradSolver.h ------------------------------*/
/*--------------------------------------------------------------------------*/

/** @file
 *Heather file for the SubGradSolver class which implements the CDASolver 
 *interface for NonDifferentiable Optimization Solvers, using a unified 
 *subgradient-type algorithm as described in:
 *
 *  A. Frangioni, E. Gorgone, B. Gendron.
 *  "On the Computational Efficiency of Subgradient Methods: a Case Study
 *   with Lagrangian Bounds"
 *  Mathematical Programming Computation 9(4), 573-604, 2017
 *
 * available at
 *
 * \link
 *  
 http://pages.di.unipi.it/frangio/abstracts.html#MPC16
 * \endlink
 *
  *
 * \version 1.0
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
 * The class requires that the function to be minimized be available under
 * the form of a C05function object, as described in C05function.h.
 *
 * The class is parametric over the type of both stepsize and direction used:
 * it just relies over the objects of classes Stepsize and Deflection.
 *
 * Copyright &copy 2001 - 2021 by Antonio Frangioni, Enrico Gorgone, Laura Fiorella
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
#include <ctime>
#include <vector>
#include <math.h>

#include "CDASolver.h"

#include "C05Function.h"
#include "LinearFunction.h"

#include "Block.h"
#include "ColVariable.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"
                        
/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{


/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SubGradSolver_CLASSES Classes in SubGradSolver.h

 *  @{ */
 class Deflection;
 class Stepsize;	//forward declaration

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SubGradSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/**
 *SubGradSolver class implements the CDASolver 
 *interface for NonDifferentiable Optimization Solvers, using a unified 
 *subgradient-type algorithm as described in:
 *
 *  A. Frangioni, E. Gorgone, B. Gendron.
 *  "On the Computational Efficiency of Subgradient Methods: a Case Study
 *   with Lagrangian Bounds"
 *  Mathematical Programming Computation 9(4), 573-604, 2017
 *
 * available at
 *
 * \link
 *  
 http://pages.di.unipi.it/frangio/abstracts.html#MPC16
 * \endlink
 *
  This is in fact a subgradient method (SM) based on abstract rules for
 * the computation of both the <em>stepsize</em> \f$\nu_i\f$ and the
 * <em>direction</em> \f$ d_i\f$. The algorithm in employs the simple
 * recurrence formula:
 *  \f[
 *    \breve{\lambda}_{i+1} \gets \bar{\lambda}_i - \nu_i d_i
 *    \quad , \quad
 *    \lambda_{i+1} \gets {\rm P}_{\Lambda}( \, \breve{\lambda}_{i+1} \, )
 *  \f]
 * where \f${\rm P}\f$ denotes the <em>orthogonal projection</em> on
 * \f$\bar{\Lambda}\f$. The point \f$\bar{\lambda}_i\f$ is not necessarily
 * the current iterate. For instance, it could be required that
 * \f$\bar{\lambda}_{i}\f$ remains unchanged if an ascent direction occurs.
 * It recalls somehow the <em>stability center</em> of the bundle methods.
 * The class relies on the objects of the friend classes <tt>Stepsize</tt>
 * and <tt>Deflection</tt>, which have to return, respectively, the
 * stepsize \f$\nu_i\f$ and the deflection coefficient \f$\alpha_i\f$. The
 * latter scalar number defines in turn the direction \f$ d_i\f$, i.e..
 * \f$ d_i = \alpha_i g_i + (1-\alpha_i)d_{i-1} \f$.
 * These abstract classes allow us to derive several variants of the method.

 * For the sake of simplicity, the original SM characterized by
 * \f$\alpha_i = 1\f$ is performed setting the pointer to the object
 * Deflection to nullptr.
 * The class also includes the <em>incremental</em> variant for when the
 * function to be maximized is composed by a sum of different functions. */
 
 
 class SubGradSolver : public CDASolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 * "Import" basic types from Function and C05Function.
 *
 *  @{ */
 
 using Index = Function::Index;
 using c_Index = Function::c_Index;

 using VarValue = Function::FunctionValue;
 using c_VarValue = Function::c_FunctionValue;

 using Vec_VarValue = Function::Vec_FunctionValue;

 using LinearCombination = C05Function::LinearCombination;
 using c_LinearCombination = C05Function::c_LinearCombination;
 
 
 
friend class Stepsize; 	// These classes are defined as friend classes, in
friend class Deflection;   // order to allow them to use directly SubgGrad's
                         // data structures [see below].
/* Note: a few methods of Stepsize and Deflection are implemented at the end of
   SubGradSolver.cpp. These are the methods that allow to read protected data of
   SubGradSolver, which may be useful to the objects. The methods are implemented
   here for the base classes, so that any derived class can use them to
   access to these data. This is crucial, as derived classes from friend
   classes are not friend, and therefore they cannot access data of SubGradSolver
   unless this capability is explicitly provided by the base classes, who are
   friends. */
 
enum int_par_type_SGSlv {
	
		/** Public enum describing the different types of algorithmic parameters
		of "int" type that any Solver should reasonably have. The value 
		intLastParCDAS is provided so that the list can be easily further 
		extended by derived classes.*/
	
	intSGPar1=CDASolver::intLastParCDAS, ///< projection-strategy parameters

	intSGPar3, ///< scheme: stepsize(deflection)-restricted

	intSGPar4,        ///< control if \f$ \hat{\lambda}_i \f$ is kept

	intSGPar5, ///< seed
								
	intRstAlgPar,	/// parameter used to choose which data of the problem you have to reset
	
	intLastParSG

		};
		
enum dbl_par_type_SGSlv {

	/** public enum for the double algorithmic parameters. Public enum describing 
	the different types of algorithmic parameters of "double" type that any Solver 
	should reasonably have. The value dblLastParCDAS is provided so that the list can
	be easily further extended by derived classes.*/

	dblSGPar2=CDASolver::dblLastParCDAS, /// incremental factor:

	dbltStar, ///[1e+2]   estimate of the longest step that can be performed

	dblEpsLin, ///[1e-6]  *relative* precision required  (doublerelacc) 

	dblEInit, /// [1e-2] 
	
	dblEFnalt, ///[1e-6] 
	
	dblEDcrs, /// [.95]  
	
	dblEStps,  /// [0]    

	dblLastParSG

		}; 
		
enum string_par_type_SGSlv {
	Stepsize_name_par=CDASolver::strLastParCDAS,
	Deflection_name_par,
	StringLastParSG
	};
/*@} -----------------------------------------------------------------------*/
/*----------------- CONSTRUCTING AND DESTRUCTING SubGradSolver --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing SubGradSolver
 *  @{ */

 /** constructor: ensure every field is initialized
  Void constructor: does nothing special */
 
 SubGradSolver( void ) : CDASolver()
 {
   Result=kUnEval;
   NumVar=0;
   NrFi=0;
   SCalls=0;
   ParIter=0;
   f_lf=nullptr;
   stepsize=nullptr;
   deflection=nullptr;
   LHasChgd=true;
   LHasProj=false;
   KpBstL=true;
   TrueLB=false;
   DoSS=true;
   GiEvaltns=0;
   CSSCntr=0;
   CNSCntr=0;
   HasLin=false;
   IsConvex=true;
   IsConcave=true;
   NrmDir=0;
   NrmGi=0;
   CSmallStep=0;
   NItIncr=0;
   FiLambda=Inf<VarValue>();
   FiHat=Inf<VarValue>();
   FiBar=Inf<VarValue>();
   FiBest=Inf<VarValue>();
   Lowerbound=-Inf<VarValue>();
   alpha=0;
   step=0;
   Sigma=0;
   Epsilon=0;
   SigmaHat=0;
   dGk=0;
   dM1Gk=0;
   dM1GkDone=false;
   ZeroComp=false;
   InnIter=false;
   DoSS=false;
   EmptySet=false;
   gi=nullptr;
    // ensure all parameters are properly given their default value
   MaxIter = CDASolver::get_dflt_int_par( intMaxIter );
   EverykIt = CDASolver::get_dflt_int_par( intEverykIt );
   LogVerb = CDASolver::get_dflt_int_par( intLogVerb );
   SGPar1=dflt_int_par[intSGPar1-intLastParCDAS];
   SGPar3=dflt_int_par[intSGPar3-intLastParCDAS];
   SGPar4=dflt_int_par[intSGPar4-intLastParCDAS];
   SGPar5=dflt_int_par[intSGPar5-intLastParCDAS];
   RstAlgPar = dflt_int_par[intRstAlgPar-intLastParCDAS];
 
   MaxTime = CDASolver::get_dflt_dbl_par( dblMaxTime );
   EveryTTm = CDASolver::get_dflt_dbl_par( dblEveryTTm );
   SGPar2=dflt_dbl_par[ dblSGPar2-dblLastParCDAS];
   tStar=dflt_dbl_par[ dbltStar - dblLastParCDAS ];
   EpsLin=CDASolver::get_dflt_dbl_par(dblRelAcc);
   EInit=dflt_dbl_par[ dblEInit - dblLastParCDAS ];
   EFnal=dflt_dbl_par[ dblEFnalt - dblLastParCDAS ];
   EDcrs=dflt_dbl_par[ dblEDcrs - dblLastParCDAS ];
   EStps=dflt_dbl_par[ dblEStps - dblLastParCDAS ];
   
   Stepsize_name=dflt_str_par[Stepsize_name_par-strLastParCDAS];
   Deflection_name=dflt_str_par[Deflection_name_par-strLastParCDAS];
 
   v_events.resize( max_event_number() );
 }
 
 /*--------------------------------------------------------------------------*/
 /// destructor: cleanly detaches the SubGradSolver from the Block

 virtual ~SubGradSolver() { set_Block( nullptr ); }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 *  @{ */
/// set the (pointer to the) Block that the Solver has to solve
void set_Block( Block * block ) override; 

/*--------------------------------------------------------------------------*/

/// set the "int" parameters of SubGradSolver
 /** Set the "int" parameters specific of SubGradSolver, together with the
  * paramaters of CDASolver that SubGradSolver actually "listens to":
  *
  * - intMaxIter [Inf<int>]: maximum iterations for the next call to solve()
  *
  * - intEverykIt [0]: after how many iteration call the eEverykIteration
  *                    events
  *
  * - intLogVerb [0]: "verbosity" of the SubGradSolver log
  *                   0 = no log
  *                   1 = only final state of the call and errors
  *                   2 = detailed step-by-step log
  *                   3 = as 2 + print every linearization added/removed
  * 
  * - intSGPar1 [0]: The direction \f$ d_i\f$ is assumed to be such a
  * 			         convex combination of the subgradient \f$ g_i\f$ and
  *				      the direction \f$ d_{i-1}\f$:
  *				      \f[
  * 			          d_i = \alpha_i g_i + ( 1 - \alpha_i ) d_{i-1}
  *				      \f]
  *				      SGPar1 selects among \f$\{d_i,d_{i-1},g_i\} \f$ the
  *					   vector(s) to be projected over the tangent cone at
  *					   \f$\bar{\lambda}_i\f$. The fields are coded bit-wise  
  *					   in the following way:
  *
  *					- SGPar1 [0]: 
  *					  if 1 the subgradient \f$ g_i\f$ is projected, 0 otherwise;
  *
  *					- SGPar1 [1]: 
  *					  if 1 the direction \f$ d_{i-1}\f$ is projected, 0 otherwise;
  *
  *					- SGPar [2]: 
  *					  if 1 the \f$ d_{i}\f$ is projected, 0 otherwise;
  *	
  * - intSGPar3 [1]: The convergence scheme. This field is coded bit-wise 
  *				      in the following way:
  *
  *					- SGPar3[0]:  1 if the safe rule is used, 0 otherwise
  *

  *					- SGPar3[1]:  1 for the stepsize-restricted scheme,
  *									0 for the deflection-restricted scheme.
  *
  * - intSGPar4 [0]: enables the use of
  *		        \f[
  *		           \hat{\lambda}_{i+1}=\alpha_{i+1}\lambda_i+(1 -\alpha_{i+1})\hat{\lambda}_i
  *		        \f]
  *		            which could have a certain influence on the stopping test 
  *					   [see IsOptimal()]. 
  *
  * - intSGPar5 [0]: The seed value used in a call to srand. The components are re-shuffled 
  *				 in the incrmental variant, and a random number generator is used.
  *
  * 
  * - intRstAlgPar [0]: Decide which dataof the problem has to be resetted. It is coded 
  *					bit-wise in the following way:
  *					
  * 				- bit 0: if 0, all the algorithmic parameters are reset to the default
  *   				  values read by the stream/set by SetPar(), while if 1 they are left
  *   				  untouched;
  *
  * 				- bit 1: if 0 the current point is reset to the all-0 vector, while if
  *   				  1 it is left untouched;
  *
  * 				- bit 2: if 0, all the subgradients are removed, except the constant
  *   				  (sub)gradient of the linear 0-th component, while if 1 the subgrads
  *   				  are left there;
  * 
  * 				- bit 3: if 0 the value of Fi() in the current point is reset to INF
  *   				  (i.e., unknown), while if 1 it is left untouched; note that resetting
  *   				  the current point [see bit 1] has this as a side-effect, regardless to
  *   				  the value of bit 3. 
  *
  **/

void set_par( const idx_type par , const int value ) override;
/*--------------------------------------------------------------------------*/
/// set the "double" parameters of SubGradSolver
 /** Set the "double" parameters specific of SubGradSolver, together with the
  * parameters of CDASolver that SubGradSolver actually "listens to":
  *
  * - dblMaxTime [Inf<double>()]: maximum CPU time for the next call to
  *                               compute(), in seconds
  *
  * - dblEveryTTm [0]: periodicity of eEveryTTime events
  *
  * - dblSGPar2 [0]: To manage the incremental iterations. A sequence
  *				      of incremental (or inner) iterations NItIncr performed
  *	  			      along single-component subgradients could occur before a
  *		            full (or normal, or outer) iteration. The number NItIncr is 
  *                  obtained as NItIncr = ceil( <Number of components of 
  *                  \f$ f\f$> * SGPar2 ), where the number of components of \f$ f\f$ 
  *                  includes the 0-th component, and then not necessarily coincides 
  *                  with NrFi. Finally, the incremental variant only works with no 
  *                  deflection object, i.e setting deflection = NULL, because we do
  *                  not manage with the sum of subgradients of different components.
  *
  * - dbltStar [1e+2]: estimate of the longest step that can be performed optimality
  *				 		   parameter related to subgradient scaling.
  *  
  * - dblEpsLin [1e-6]: Optimality related parameters. 
  *					Proving that some point Lambda is optimal for a NonDifferentiable
  *					Optimization problem involves finding an all-0  subgradient	of the
  *					function at Lambda. If an all-0 vector is found in the 
  *					epsilon-subdifferential of Lambda, then the point is epsilon-optimal.
  *					Note that if the minimization problem is subject to constraints, i.e.,
  *					Fi() has to be minimized only on the points Lambda \in L, the latter
  *					being a convex set, then the above is referred to a subgradient of the
  *					"actual function" ( Fi + I_L )( Lambda ), where I_L is the indicator
  *					function of L (evaluating to 0 inside L and to +INF otherwise). In other
  *					words, one has to show that there exists a( epsilon-)subgradient of
  *					Fi() at Lambda that, *after projection on the frontier of L*, is all-0.
  *					A general stopping condition requires that, if EpsLin is the *relative*
  *					precision required, a solver can stop if it finds an epsilon-subgradient
  *					g at Lambda such that
  *								tStar * || g || + epsilon <= EpsLin * | MaxFi |
  *					where MaxFi is an estimate of the optimal solution value of the NDO
  *					problem, tStar is an estimate of the longest step that can be performed
  *					and || || is a norm-like function. tStar is related to the "scaling" of
  *					Fi(), and it can be seen as an estimate of the actual decrease that can be
  *					obtained by moving of an unitary step in the direction of any subgradient.
  *					Alternatively, the above condition can be seen as a weaker form of
  *						epsilon <= EpsLin*| MaxFi |/2|| g ||<=EpsLin*| MaxFi |/(2*tStar)
  *					which says that the solver stops when both epsilon and || g || are
  *					"small", with tStar dictating what "small" means for || g ||. Note that
  *					each derived class can use different norm-like functions to evaluate one
  *					or the other of the above conditions. 
  *
  * - dblEInit [1e-2] 
  *	
  *	- dblEFnalt [1e-6] 
  *	
  *	- dblEDcrs [0.95]  
  *	
  *	- dblEStps [0]: The evaluation of function to be minimized may be a costly task: in many  
  *					cases, it requires the solution of a - possibly hard - optimization problem. 
  *					Often, time can be saved if the function is only approximately computed at 
  *					the beginning of the optimization process; of course, the computation should 
  *					become more and more "exact" as the optimization proceeds. The relative 
  *					precision required is initially set to EInit, and decreased down to EFnal 
  *					by multiplying it by EDcrs every EStps iterations. If the algorithm allows it,
  *					EStps can be set to 0 meaning that the precision is decreased only if necessary,
  *					i.e., when it is impossible to proceed otherwise (the algorithm must have some 
  *					way to detect this). The precision is kept fixed if EInit == EFnal. A solver
  *				    which is *not* capable of handling approximate computation of the function
  *				    can ignore the values of these parameters.
  *
*/

void set_par( const idx_type par , const double value ) override;
 /*--------------------------------------------------------------------------*/
 /// set the string parameters
 void set_par (idx_type par , std::string && value ) override;

 /*--------------------------------------------------------------------------*/
 /// set the ostream for the SubGradSolver log

 void set_log( std::ostream * log_stream = nullptr ) override;

/*@} ----------------------------------------------------------------------*/
/*----------------- METHODS FOR ACCESSING THE DATA OF THE Block ------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the data of the Block
 *
 * These methods provide convenient shortcuts for directly asking to the
 * SubGradSolver some relevant data about the Block it is solving.
 *
 *  @{ */
/// returns the number of "components", i.e., C05Function in the objective
 /** Returns the total number of "components", i.e., the C05Function whose
  * sum (possibly together with one LinearFunction) makes up the objective of
  * the Block that the SubGradSolver is solving. This method should not be
  * called if set_Block() has not been called, or has last been called with
  * nullptr argument. */

 Index n_components( void ) const { return( NrFi ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns a pointer to the i-th C05Function in the objective
 /** Returns a pointer to the i-th C05Function, i.e., the i-th term of the
  * sum of C05Function that (possibly together with one LinearFunction) 
  * makes up the objective of the Block that the SubGradSolver is solving. 
  * It must ve 0 <= \p i <= n_components(). All returned pointers are not
  * nullptr provided that set_Block() has last been called with not nullptr
  * argument (otherwise this method should not be called). */

 C05Function * component( Index i ) const {
  #ifndef NDEBUG
   if( i >= NrFi )
    throw( std::invalid_argument( "wrong component number" ) );
  #endif
  return( v_c05f[ i ] );
  }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns a pointer to the LinearFunction in the objective
 /** Returns a pointer to the single LinearFunction that is summed with the
  * C05Function in the objective, if any. If the method returns nullptr,
  * there is no such LinearFunction (it is constantly 0, i.e., all its
  * coefficients are 0). This method should not be called if set_Block() has
  * not been called, or has last been called with nullptr argument. */

 LinearFunction * l_component( void ) const { return( f_lf ); }
  
  /**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR EVENTS HANDLING -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Set event handlers
 *
 *  SubGradSolver heeds to all three "basic" types of events:
 *
 * - eBeforeTermination, called just before optimality stop (but not all
 *   other kinds of stop), with return action eForceContinue forcing one
 *   new iteration (master problem solution) to be performed;
 *
 * - eEverykIteration, called every intEverykIt iterations (if intEverykIt
 *   != 0), with possible return actions eStopOK and eStopError;
 *
 * - eEveryTTime, called every dblEveryTTm seconds (if dblEveryTTm != 0),
 *   with possible return actions eStopOK and eStopError;
 *
 * Of course, events have been set with set_event_handler() for them to be
 * called.
 *  @{ */

 EventID max_event_number( void ) const override { return( 3 ); }

 /*@} 
----------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the NDO encoded by the current Block
 *  @{ */

 /// (try to) solve 

 int compute( bool changedvars = true ) override;

/**< Tries to minimize the function. It implements the subgradient algorithm
 * exploiting the protected methods FormD(), SaveDir(), FormLambda1(),
 * FiAndGi(), and GotoLambda1().
 *
 * Returns   if
 *
 * - kOK       optimization has been succesful: a solution that is "optimal"
 *             (w.r.t. the current parameters settings) has been found;
 *
 * - kUnbunded   there has been an error: Fi() has returned
 *             -INF, or the function is unbounded below: the latter case can
 *             be detected only if a lower bound on the min. value of
 *             \f$f\f$ is available;
 *
 * - kUnfesible   the polyhedral set defined by the constraints is empty: in this
 *             case, the primal optimal solution is an unbounded *extreme ray*
 *             for the dual problem;
 *
 * - kStopped  compute() has been stopped, either by Function::compute()
 *             or because the stepsize has been "too small" during 100
 *             consecutive outer iterations [or 100 * NrFi consecutive inner
 *             iterations];
 *
 * - kStopIter  the max. number of iterations has been exhausted;
 *
 * - kStopTime  the max. running time has been reached;
 *
 * - kError    There have been some problem in the Block that require to
 *             stop the optimization.
 *
 * As for kStopped, "too small" means that \f$ \nu_i \leq 1e-8 * t^*\f$, where
 * \f$ t^* \f$ is the optimality related parameter scaling Fi(). There is no
 * reason, in principle, why we couldn't replace \f$1e-8\f$ by a parameter,
 * but in order to make the test easier this parameter has been fixed to
 * \f$1e-8\f$. We also decided to replace by 100 the parameter saying how many
 * outer iterations of consecutive small stepsizes are sufficient to stop the
 * algorithm.
 *
 * Note that, whatever the exit condition be, the current point is always
 * available by calling ReadSol(), and its Fi() value by calling ReadFiVal().
 */
/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

/// get the number of Variable
/** Variable cannot be changed. This means that is used the default
 *  implementation of GetMaxNumVar(). The maximum number of variables is
 *  equal to the current number of variable*/

   Index GetNumVar( void );

/*--------------------------------------------------------------------------*/

   Index GetNrFi( void );

/*--------------------------------------------------------------------------*/
/** Returns the pointer to the function */

   std::vector< C05Function * > GetFun( void );

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
@{ */

   VarValue ReadBestFiVal( void );

/**< Returns the best \f$f\f$ value found so far. Independently from which
 * "component" of Fi() is chosen, it returns the full function. */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  Vec_VarValue ReadBestSol( void );

   Vec_VarValue ReadSol( void );

/**< Returns a vector of the <em>stability center</em>
 * \f$\bar{\lambda}_i\f$. If compute() has returned a kOK and the tStar has
 * been properly set, the point returned by ReadSol() is \f$\epsilon\f$-optimal. */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   VarValue ReadFiVal( void );

/**< Independently from which "component" of Fi() is chosen, it returns the

 * full function Fi at the <em> stability center </em> \f$\bar{\lambda}_i\f$.
 */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  VarValue ReadHatFiVal( void );

/**< Returns \f$\hat{f}\f$, if \f$ \hat{\lambda}_i \f$ is kept in memory.
 * Otherwise it returns Inf<HpNum>().   */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /*--------------------------------------------------------------------------*/  
 /** Returns the current number of iterations. Note that this counter is *not*
    necessarily reset each time compute() is called, as multiple calls may
    actually be necessary to "solve the given problem". Thus, ideally this
    counter is only reset when "it is clear that an entirely new problem is
    being solved", such as after a call to ReSetAlg() or when the Block is changed. */

   Index NrIter( void ) const { return( ParIter ); }
   
/*--------------------------------------------------------------------------*/  

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the elapsed CPU time since the start of the last compute()

 double elapsed( void ) const {
  return( ( std::clock() - c_start) / CLOCKS_PER_SEC );
  }
  

bool IsOptimal( VarValue eps = 0 );

/**< Returns true if the solution \f$\bar{\lambda}_i\f$ is
 * \f$\epsilon\f$-optimal (relative), being \f$\epsilon \f$ set to EpsLin.
 * The parameter SGPar4 controls the linearization error to be used in 
 * the stopping condition:
 *
 * -# SGPar4 = true:
 *    \f[
 *     t^* \|d_i\| + \max\{ \hat{\epsilon}_i , \epsilon_i \} <= 
 *     \epsilon * max( 1 , |f_i^{rec}| )
 *    \f]
 *
 * -# SGPar4 = false: The criteria is just
 *    \f[
 *     t^* \|d_i\| + \epsilon_i <= \epsilon * max( 1 , |f_i^{rec}| )
 *    \f]
 *
 * where \f$ f^{rec}_i = \min \{ \, f_l \,:\, l = 1, \ldots, i \, \}\f$, i.e. 
 * the <em> record value </em> on the optimum $f_*$ */

/*--------------------------------------------------------------------------*/

/**@} ----------------------------------------------------------------------*/
/*************************************************/
/** @name Accessing the found solutions (if any)
 *  @{ */

  VarValue get_lb( void ) override 
  {
   if(IsConcave)
	   return (FiBest ? -FiBest : -Inf<VarValue>());
	else{
		if(IsConvex)
			return(TrueLB ? Lowerbound :-Inf<VarValue>());
		else
			throw(std::invalid_argument("The function is nieghter concave or convex"));
	}
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 VarValue get_ub( void ) override
 {
   if(IsConvex)
	   return (FiBest ? FiBest : Inf<VarValue>());
	else{
		if(IsConcave)
			return(TrueLB ? -Lowerbound :Inf<VarValue>());
		else
			throw(std::invalid_argument("The function is nieghter concave or convex"));
	}
  }

/*--------------------------------------------------------------------------*/

 /// write the "current" solution

 void get_var_solution( Configuration *solc = nullptr ) override{
  /*if( ( MaxSol > 1 ) && ( UpFiBest < UpRifFi[ NrFi ] ) ) {
   for( Index i = 0 ; i < NumVar ; i++ )
    LamVcblr[ i ]->set_value( LmbdBst[ i ] );
   }
  else {
   for( Index i = 0 ; i < NumVar ; i++ )
    LamVcblr[ i ]->set_value( Lambda[ i ] );*/
    throw("Not implemented yet");
   }


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// write the "current" dual solution

 void get_dual_solution( Configuration *solc = nullptr ) override{
 throw("Not implemented yet");
 }

/*--------------------------------------------------------------------------*/

 /*bool new_var_solution( void ) override
 {
  if( ( MaxSol > 1 ) && ( UpFiBest < UpRifFi[ NrFi ] ) ) {
   // dirty trick: pretend that LmbdBst is not better than Lambda by
   // re-defining UpFiBest, so that next time Lambda is given
   // anyway, the value of UpFiBest is no longer used after that the

   // corresponding solution has been given and a new one if "generated"

   UpFiBest = UpRifFi[ NrFi ];

   return( true );
   }
  else
   return( false );
  }
  */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 //bool new_dual_solution( void )  override { return( false ); }

/*--------------------------------------------------------------------------*/


 /*@} -----------------------------------------------------------------------*/
 /*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
	/** @name Handling the parameters of the Solver
	     *  @{ */

 virtual idx_type get_num_int_par( void ) const override{
	      return( idx_type( intLastParSG ) );
	      }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual idx_type get_num_dbl_par( void ) const override{
	      return( idx_type( dblLastParSG ) );
	      }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual idx_type dbl_par_str2idx( const std::string & name ) const override{
	    	 // these may be many enough as to warrant using a map
	    	 const auto it = dbl_pars_map.find( name );
	    	 if( it != dbl_pars_map.end() )
	    	  return( it->second );
	    	 else
			   return( Solver::dbl_par_str2idx( name ) );
	}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	     
virtual idx_type int_par_str2idx( const std::string & name ) const override{
	    	  // these may be many enough as to warrant using a map
	    	  const auto it = int_pars_map.find( name );
	    	  if( it != int_pars_map.end() )
	    	   return( it->second );
	    	  else

	       return( Solver::int_par_str2idx( name ) );
 }
 
 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	     
virtual idx_type str_par_str2idx( const std::string & name ) const override{
	    	  // these may be many enough as to warrant using a map
	    	  const auto it = str_pars_map.find( name );
	    	  if( it != str_pars_map.end() )
	    	   return( it->second );
	    	  else

	       return( Solver::str_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

virtual const std::string & dbl_par_idx2str( const idx_type idx )const override{

	    	  	        if( ( idx >= dblSGPar2 ) && ( idx < dblLastParSG ) )
	    	  	         return( dbl_pars_str[ idx - dblSGPar2 ] );
	    	  	        else
	    	  	         return( Solver::dbl_par_idx2str( idx ) );
}
/*--------------------------------------------------------------------------*/

 virtual int get_int_par( const idx_type par ) ;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 virtual double get_dbl_par( const idx_type par ) ;
 
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 //virtual double get_str_par( const idx_type par ) ;

 /*--------------------------------------------------------------------------*/
void ReSetAlg( unsigned char RstLvl=0 );

 /* Resets the internal state of the subgrad algorithm. Since several
    different things can be reset independently, RstLvl is coded bit-wise:

    - bit 0: if 0, all the algorithmic parameters are reset to the default
      values read by the stream/set by SetPar(), while if 1 they are left
      untouched;

    - bit 1: if 0 the current point is reset to the all-0 vector, while if
      1 it is left untouched;

    - bit 2: if 0, the current point is reset to the value currently in the
      active Variable of the C05Function, while if 1 it is left untouched. */
	  
/*--------------------------------------------------------------------------*/
virtual const std::string & int_par_idx2str( const idx_type idx )const override{

	        if( ( idx >= intSGPar1 ) && ( idx < intLastParSG ) )
	         return( int_pars_str[ idx - intSGPar1 ] );
	        else
	         return( Solver::int_par_idx2str( idx ) );
}

/*--------------------------------------------------------------------------*/
	        /// get the default value of an int parameter
	        /** Get the default value of the int parameter with given index. The method
	         * is given a void implementation (returning 0), rather than being pure
	         * virtual, so that derived classes not having any int parameter do not
	         * have to bother with implementing it. */

virtual int get_dflt_int_par( const idx_type par ) const override{
	        	
            if( ( par >= intSGPar1 ) && ( par < intLastParSG ) )
	        	  return( dflt_int_par[ par - intSGPar1 ] );
	        	else
	        	  return( CDASolver::get_dflt_int_par( par ) );
}

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	        /// get the default value of a double parameter
	        /** Get the default value of the double parameter with given index. The
	         * method is given a void implementation (returning 0), rather than being
	         * pure virtual, so that derived classes not having any double parameter do
	         * not have to bother with implementing it. */

virtual double get_dflt_dbl_par( const idx_type par ) const override{
	        	
            if( ( par >= dblSGPar2 ) && ( par < dblLastParSG ) )
	        	  return( dflt_dbl_par[ par - dblSGPar2 ] );
	        	else
	        	  return( CDASolver::get_dflt_dbl_par( par ) );
}

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	        /// get the default value of a string parameter
	        /** Get the default value of the string parameter with given index. The
	         * method is given a void implementation (returning the empty string),
	         * rather than being pure virtual, so that derived classes not having any
	         * string parameter do not have to bother with implementing it. */

/*virtual std::string & get_dflt_str_par( const idx_type par ) const override {
	         if( (par>=Stepsize_name_par) && (par< StringLastParSG) )
	         	return(dflt_str_par[par-Stepsize_name_par]);
	         else
	         	return( CDASolver::get_dflt_str_par(par));
}*/
 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
/*@} -----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
void SetLambda(Vec_VarValue tLambda);

void FormD( void );
		  /**< The method is used within compute(), and its job is to compute the
	   	  * direction \f$d_i\f$ that appears in the formula of \f$\lambda_{i+1}\f$.
	   	  * The direction is actually saved after the variables generation because i)
	   	  * the subgradient keeps in memory just the current direction, and ii) by
	   	  * generating/removing variables the direction may quickly come to be
	   	  * deteriorated. When the variables generation is ended , SaveDir() must
	   	  * be called in order to update the direction.
	   	  * In addition, the deflection coefficient is computed inside FormD(). As for
	   	  * the stepsize, the computation is performed within FormD() only if the
	   	  * scheme is <em>deflection-restricted</em>. */

/*--------------------------------------------------------------------------*/
 void SaveDir( void );

	   	 /**< The method is combined with FormD() to carry out the
	   	  * direction computation. Moreover, the stepsize is determined into SaveDir()
	   	  * when the adopted scheme is <em>stepsize-restricted</em>. */

/*--------------------------------------------------------------------------*/

 void FiAndGi( const Variable::Index wFi = Inf<Variable::Index>() );

	     /**< Evaluates the function at the new point \f$\lambda_{i+1}\f$, i.e.,
	      * \f$f(\lambda_{i+1})\f$, and it either computes a subgradient
	      * \f$g_{i+1} \in \partial f(\lambda_{i+1})\f$, or, if the point is
	      * infeasible, a constraint. */

/*--------------------------------------------------------------------------*/

void FormLambda1( void ); 

	   /**< After a (succesfull) call to FormD(), sets the new (unprojected)
	    * tentative point \f$\breve{\lambda}_{i+1}\f$  as
	    * \f[
	    *  \breve{\lambda}_{i+1} = \lambda_i - \nu_i d_i \,.
	    * \f]
	    * Remark that the point \f$\breve{\lambda}_{i+1}\f$ must be projected before
	    * calling FiandGi() [see above], i.e., \f$ \lambda_{i+1} = {\rm P}_{\Lambda}
	    * ( \breve{\lambda}_{i+1} ) \f$. */

/*--------------------------------------------------------------------------*/

void GotoLambda1( void );
		/**< Move the current point to \f$\lambda_{i+1}\f$.*/

/*--------------------------------------------------------------------------*/

Index GetGi(VarValue *sub, Index  * sgb, Index strt, Index stp); 
      /**< Compute the subgradient if the function has diagonal linearization*/

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 int SGPar1;   ///< projection-strategy parameter
 
 int SGPar3;   ///< scheme: stepsize(deflection)-restricted
 
 int SGPar4;   ///< control if \f$ \hat{\lambda}_i \f$ is kept
 
 int SGPar5;   ///< seed
 
 double SGPar2;   ///< incremental factor
 
 double tStar;    ///< estimate of the longest step that can be performed
 
 double EDcrs;    ///< *relative* precision required  (dblRealAcc) 
 
 double EFnal;
 
 double EInit;
 
 double EStps;
 
 double EpsLin;

 int RstAlgPar;     ///< reset parameter, bit-wise coded

 int LogVerb;       ///< "verbosity" of the log

 Index MaxIter;     ///< maximum number of iterations
 
 double MaxTime;    ///< maximum time (in seconds) for each call to Solve()

 std::string Stepsize_name;    ///< name of the Stepsize method

 std::string Deflection_name;    ///<name of the deflection method ("null" if no deflection is used)
 // generic parameter - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int Result;        ///< result of the latest call to compute()

Index NumVar;      ///< (current) number of variables
 
Index NrFi;        ///< number of components of Fi()

Index SCalls;      ///< number of calls to compute() (the current included)

Index ParIter;     ///< number of iterations in this call to compute()

Index GiEvaltns;   ///< total number of Gi() calls

Stepsize * stepsize; ///<pointer to the Stepsize class
 
Deflection *deflection;   ///< pointer to the Deflection class

Vec_VarValue Lambda;   ///< the current point

Vec_VarValue LambdaHat;  ///< the point \f$ \hat{\lambda}_i \f$

Vec_VarValue LambdaBar;     ///< the stability center \f$\bar{\lambda}_i\f$

Vec_VarValue LmbdBst;  ///< the best point found so far
 
VarValue FiLambda;      ///< full function value at \a Lambda

VarValue FiHat;         ///< full function value at \a LambdaHat

VarValue FiBar;         ///< full function value at \a LambdaBar

bool HasLin; ///<true if the function has diagonal linearization>
 
bool LHasChgd;       /**< true if Lambda has changed since the latest call
		       * to FiAndGi(): allows repeated calls in the same
		       * Lambda, e.g. with increasing precision */
			   
bool LHasProj;       /**< true if \a Lambda  has projected in the current
                         * iteration */

		 
bool KpBstL; ///<true if LmbdBst has to be kept> 

VarValue Lowerbound;    ///< Lower Bound over the full function \f$f\f$

bool TrueLB;         /**< true if LowerBound is a "true" lower bound rather
                         * than just the "minus infinity" */

Vec_VarValue Gi;            ///< Gi[ wFi ]( Lambda ), the subgradient vector
VarValue * gi;
Vec_VarValue dir;           ///< the direction \f$ d_i \f$

VarValue alpha;         ///< the deflection parameter \f$ \alpha_i \f$
   
VarValue step;          ///< the stepsize \f$ \nu_i \f$

VarValue Sigma;         /**< the linearization error \f$ \sigma_i\f$ of
                         * \f$ g_i\f$ at \f$\bar{\lambda}_{i}\f$ */

VarValue Epsilon;       /**< the linearization error \f$\epsilon_i\f$ of
                         * \f$d_i\f$ at \f$\bar{\lambda}_i\f$ */

VarValue SigmaHat;      /**< the linearization error \f$ \hat{\alpha}_i \f$ of
                         * \f$g_i\f$ with respect to \f$\hat{\lambda}_i\f$ */
VarValue HatEpsilon;    /**< the linearization error \f$ \hat{\epsilon}_i \f$
			 * of \f$d_i\f$ with respect to \f$\hat{\lambda}_i\f$
			 */
std::vector<Index> SGBase;	///< the set of indices of Gi[ wFi ]( Lambda )

VarValue NrmGi;         ///< the (squared) subgradient's norm

VarValue NrmDir;        ///< the (squared) direction's norm

VarValue dGk;           ///< scalar product \f$ < d_i , g_i > \f$
   
VarValue dM1Gk;         ///< scalar product \f$ < d_{i-1} , g_i > \f$

bool dM1GkDone;      /**< true if the scalar product
			   \f$d_{i-1}^{\top} g_i\f$ has been computed */

bool ZeroComp;       ///< true if Fi() comes with the 0-th component

Index NItIncr;       /**< number of incremental iterations after an outer
			   iteration */
bool InnIter;        ///< if true, the current iteration is an inner one

Index CSSCntr;       ///< counter of consecutive SS

Index CNSCntr;       ///< counter of consecutive NS

Index CSmallStep;    ///< counter of consecutive \a short step

bool DoSS;           ///< SS vs NS

bool EmptySet;       ///< true, if the feasible set is empty

double EveryTTm;   ///< periodicity of eEveryTTime events

int EverykIt;      ///< periodicity of eEverykIteration events

/* -------------------------------------------------------------------------------------------------*/
 VarValue FiBest;        ///< the best value of \f$f\f$ found so far

 std::vector< C05Function * > v_c05f; /**< the vector of the components of the
                                       * sum function */
 LinearFunction * f_lf;  ///< the 0-th component of the sum function
 
 std::vector< ColVariable * > LamVcblr;  ///< map Lambda -> ColVariable
 
 bool IsConvex; /**< true if evrey component of Fi() and Fi() iteself are convex function,
 		* default is <true>. Cointais the result of Function::is_convex() >*/
 
 bool IsConcave;
 int fs; /// return status from C05Function::compute()
 
 Vec_VarValue Ubv;	///< upper bounds on the variables
 
 Vec_VarValue Lbv;	///< lower bounds on the variables

 
 // fields for events - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::clock_t c_start;   ///< starting instant of last call to compute()

 // static fields - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 const static std::vector<int> dflt_int_par;
 ///< the (static const) vector of int parameters default values

 const static std::vector<double> dflt_dbl_par;
 ///< the (static const) vector of double parameters default values
 
 const static std::vector<int> dflt_str_par;
 ///< the (static const) vector of string parameters default values

 const static std::vector< std::string > int_pars_str;
 ///< the (static const) vector of int parameters names

 const static std::vector< std::string > dbl_pars_str;
 ///< the (static const) vector of double parameters names

const static std::vector< std::string > str_pars_str;
 ///< the (static const) vector of string parameters names
 const static std::map< std::string , idx_type > int_pars_map;
  ///< the (static const) map for int parameters names

 const static std::map< std::string , idx_type > dbl_pars_map;
 ///< the (static const) map for double parameters names

 const static std::map< std::string , idx_type > str_pars_map;
 ///< the (static const) map for string parameters names

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:
 
 const int FirstintDeflectionPar=intLastParSG;
 
 const int FirstintStepsizePar=FirstintDeflectionPar+100;
 
 const int FirstdblDeflectionPar=dblLastParSG;
 
 const int FirstdblStepsizePar=FirstdblDeflectionPar+100;
 /*--------------------------------------------------------------------------*/

 void guts_of_destructor( void );
 /*--------------------------------------------------------------------------*/
  
 void Log1( void );

 void Log2( void );

 void Log3( void );

 void Log4( void );


/*--------------------------------------------------------------------------*/

//inline VarValue lbnd( c_Index i ) {
	// if lb is defined return lb[i], else return -Inf
  //  return( lb ? lb[ i ] : - Inf<Vec_VarValue>() );
 //   }
//probabilmente servono
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

//inline Vec_VarValue ubnd( Index i ) { 
// if ub is defined return ub[i], else return +Inf
//	if(IsConvex)	
//		return( lb ? lb[ i ] : Inf<Vec_VarValue>() );
//	else
//		throw("Upperbound in maximization problem(?) ");
// }

/*--------------------------------------------------------------------------*/

void EvalGi( Index wFi = Inf<Index>() );

/*--------------------------------------------------------------------------*/

bool ProjFsb( Vec_VarValue Lpt , bool & EmptySet );

void ProjTng( Vec_VarValue Gpt , Vec_VarValue Lpt , c_Index strt = 0 , Index stp = Inf<Index>() );

VarValue ScalarProduct( Vec_VarValue A, Vec_VarValue B, Index n);

VarValue Norm(Vec_VarValue A, Index n);
/*--------------------------------------------------------------------------*/

void InitStepsize( void );

void InitDeflection( void );

/*--------------------------------------------------------------------------*/

void ChgLambdaHat( void );

/*--------------------------------------------------------------------------*/

void UpdateSigma( void );

/*--------------------------------------------------------------------------*/

void UpdtLowerBound( void );

/*--------------------------------------------------------------------------*/

void AllocateUBLB( c_Index strt , c_Index num );

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
};//end Class SubGradSolver

};//end namespace SMSpp_di_unipi

#endif //SubGradSolver.h included

/*--------------------------------------------------------------------------*/
/*------------------------- End File SubGradSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
