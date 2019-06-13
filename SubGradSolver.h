/*--------------------------------------------------------------------------*/
/*-------------------------- File SUBGRADSOLVER.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Definition of the SubGradSolver class , which implements the CDASolver interface
 * for NonDifferentiable Optimization Solvers, as described in CDASolver.h,
 * using a (deflected, projected, incremental) subgradient-type algorithm.
 *
 * The user is assumed to be familiar with the algorithm: refer to
 *
 *  A. Frangioni, E. Gorgone, B. Gendron.
 *  "On the Computational Efficiency of Subgradient Methods: a Case Study
 *   with Lagrangian Bounds"
 *  Mathematical Programming Computation 9(4), 573-604, 2017
 *
 * \version 1.0
 *
 * \date 27 - 04 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Department of Informatics \n
 *         University of Pisa \n\n
 *
 * \author Enrico Gorgone \n
 *         Department of Mathematics and Informatics \n
 *         University of Cagliari \n\n
 *
 *
 * \author Simone Zanda \n
 *         Department of Mathematics and Informatics \n
 *         University of Cagliari \n\n
 *
 * The class requires that the function to be minimized be available under
 * the form of a C05function object, as described in C05function.h.
 *
 * The class is parametric over the type of both stepsize and direction used:
 * it just relies over the objects of classes Stepsize and Deflection.
 *
 * Copyright &copy 2001 - 2019 by Antonio Frangioni, Enrico Gorgone,Simone Zanda
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/


#ifndef SUBGRADSOLVER_H_
#define SUBGRADSOLVER_H_

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/





#include "CDASolver.h"
#include "Function.h"
#include "C05Function.h"
#include "LinearFunction.h"

#include "Block.h"
#include "ColVariable.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"

#include "Configuration.h"

#include "Deflection.h"
#include "Stepsize.h"
















/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SubGrad_CLASSES Classes in SubGrad.h
    @{ */

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS SubGradSolver  -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** The SubGradSolver class implements the CDASolver interface for NonDifferentiable
 * Optimization Solvers, using a unified subgradient-type algorithm as
 * described in:
 *
 *  A. Frangioni, E. Gorgone, B. Gendron.
 *  "On the Computational Efficiency of Subgradient Methods: a Case Study
 *   with Lagrangian Bounds"
 *  Mathematical Programming Computation 9(4), 573-604, 2017
 *
 * This is in fact a subgradient method (SM) based on abstract rules for
 * the computation of both the <em>stepsize</em> \f$\nu_i\f$ and the
 * <em>direction</em> \f$d_i\f$. The algorithm in employs the simple
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
 * latter scalar number defines in turn the direction \f$d_i\f$, i.e..
 * \f$ d_i = \alpha_i g_i + (1-\alpha_i)d_{i-1} \f$.
 * These abstract classes allow us to derive several variants of the method.
 * For the sake of simplicity, the original SM characterized by
 * \f$\alpha_i = 1\f$ is performed setting the pointer to the object
 * Deflection to nullptr.
 * The class also includes the <em>incremental</em> variant for when the
 * function to be maximized is composed by a sum of different functions. */



class Stepsize; 	// forward declaration of class Stepsize
class Deflection;	// forward declaration of class Deflection


class SubGradSolver : public CDASolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/



public:

	/*--------------------------------------------------------------------------*/
	/*---------------------------- PUBLIC TYPES --------------------------------*/
	/*--------------------------------------------------------------------------*/
	/** @name Public Types
	    @{ */



	//Public enum describing the "status" of the FiOracle as returned by
		  //  GetFiStatus() [see below].

		   enum FiStatus { kFiNorm = 0,
				   kFiError ,
				   kFiChgd ,
				   kFiStop ,
				   kFiCont
		                   };

		 /// Public enum "extending" sol_type to a specific case of CDASolvers

		   enum sol_type {
		     kEILoopNow = kBothInfeasible + 1 ,

		     kEIAbort,
		     };



		enum int_par_type_SG{

			 ///< first allowed new int parameter for derived classes
		    /**< Convenience value for easily allow derived classes
		* to extend the set of int algorithmic parameters. */

			intSGPar1=intLastAlgPar, ///< projection-strategy parameters

			/*-#  SGPar1 [0] The direction \f$d_i\f$ is assumed to be such a
				                    convex combination of the subgradient \f$g_i\f$ and
				                    the direction \f$d_{i-1}\f$:
				                    \f[
				                      d_i = \alpha_i g_i + ( 1 - \alpha_i ) d_{i-1}
				                    \f]
				                    SGPar1 selects among \f$\{d_i,d_{i-1},g_i\} \f$ the
				       vector(s) to be projected over the tangent cone at
				       \f$\bar{\lambda}_i\f$. The field SGPar1 is coded bit-wise, in the
				       following way:

				        - bit 0:  if 1 the subgradient \f$g_i\f$ is projected, 0 otherwise;

				        - bit 1:  if 1 the direction \f$d_{i-1}\f$ is projected, 0 otherwise;

				        - bit 2:  if 1 the \f$d_{i}\f$ is projected, 0 otherwise;

				       The setting (+7) is redundant. In fact, it is equivalent to (+3)
				       because \f$d_i\f$, being a (convex) combination of \f$g_i\f$ and
				       \f$d_{i-1}\f$, coincides with its projection.*/

			intSGPar3, ///< scheme: stepsize(deflection)-restricted

			/* -# Index SGPar3 [1] The convergence scheme. This field is coded
		                        bit-wise in the following way:

		                         - bit 0:  1 if the safe rule is used, 0 otherwise

		                         - bit 1:  1 for the stepsize-restricted scheme,
		                                   0 for the deflection-restricted scheme.*/

			intSGPar4,        ///< control if \f$ \hat{\lambda}_i \f$ is kept

			/*-# bool SGPar4 [true] SGPar4 enables the use of
		                          \f[
		                           \hat{\lambda}_{i+1} = \alpha_{i+1}\lambda_i +
		                           ( 1 - \alpha_{i+1} ) \hat{\lambda}_i
		                          \f]
		                          which could have a certain influence on the
		                          stopping test [see IsOptimal()].*/

			intSGPar5, ///< seed

			/* -# Index SGPar5 [0] The seed value used in a call to srand. The
		                        components are re-shuffled in the incrmental variant,
		                        and a random number generator is used.*/


			intLastAlgParSG,



		};// end( int_par_type_SG )


		enum dbl_par_type_SG {

			   ///< first allowed new double parameter for derived classes
			                 /**< Convenience value for easily allow derived classes
					  * to extend the set of double algorithmic parameters. */

			dblSGPar2=dblLastAlgPar, ///< incremental factor

			/* -#  SGPar2 [0] To manage the incremental iterations. A sequence
		                    of incremental (or inner) iterations NItIncr performed
		                    along single-component subgradients could occur before a
		       full (or normal, or outer) iteration. The number NItIncr is obtained
		       as NItIncr = ceil( <Number of components of \f$f\f$> * SGPar2 ),
		       where the number of components of \f$f\f$ includes the 0-th component,
		       and then not necessarily coincides with NrFi. Finally, the incremental
		       variant only works with no deflection object, i.e setting deflection
		       = NULL, because we do not manage with the sum of subgradients of
		       different components.*/

			dbltStar, //[1e+2]   Optimality related parameters. Proving that some

			dblEpsLin, //[1e-6] point Lambda is optimal for a NonDifferentiable

			/*Optimization problem involves finding an all-0 subgradient
			   of the function at Lambda. If an all-0 vector is found in
	    the epsilon-subdifferential of Lambda, then the point is epsilon-optimal.
	    Note that if the minimization problem is subject to constraints, i.e.,
	    Fi() has to be minimized only on the points Lambda \in L, the latter
	    being a convex set, then the above is referred to a subgradient of the
	    "actual function" ( Fi + I_L )( Lambda ), where I_L is the indicator
	    function of L (evaluating to 0 inside L and to +INF otherwise). In other
	    words, one has to show that there exists a( enspilon-)subgradient of
	    Fi() at Lambda that, *after projection on the frontier of L*, is all-0.
	    A general stopping condition requires that, if EpsLin is the *relative*
	    precision required, a solver can stop if it finds an epsilon-subgradient
	    g at Lambda such that
	              tStar * || g || + epsilon <= EpsLin * | MaxFi |
	    where MaxFi is an estimate of the optimal solution value of the NDO
	    problem, tStar is an estimate of the longest step that can be performed
	    and || || is a norm-like function. tStar is related to the "scaling" of
	    Fi(), and it can be seen as an estimate of the actual decrease that can be
	    obtained by moving of an unitary step in the direction of any subgradient.
	    Alternatively, the above condition can be seen as a weaker form of
	        epsilon <= EpsLin * | MaxFi | / 2
		|| g || <= EpsLin * | MaxFi | / ( 2 * tStar )
	    which says that the solver stops when both epsilon and || g || are
	    "small", with tStar dictating what "small" means for || g ||. Note that
	    each derived class can use different norm-like functions to evaluate one
	    or the other of the above conditions. Also, different estimates can be
	    used for MaxFi, although using the bast Fi-value found so far is pretty
	    common.*/

			dblEInit, // [1e-2] The evaluation of function to be minimized may be a
			dblEFnalt, //[1e-6] costly task: in many cases, it requires the solution
			dblEDcrs, // [.95]  of a - possibly hard - optimization problem. Often,
			dblEStps,  // [0]    time can be saved if the function is only
			/*approximately computed at the beginning of the optimization
			   process; of course, the computation should become more and
	    more "exact" as the optimization proceeds. The relative precision required
	    to the FiOracle [see SetPrecision() in FiOracle.h] is initially set to
	    EInit, and decreased down to EFnal by multiplying it by EDcrs every EStps
	    iterations. If the NDO algorithm allows it, EStps can be set to 0 meaning
	    that the precision is decreased only if necessary, i.e., when it is
	    impossible to proceed otherwise (the algorithm must have some way to
	    detect this). The precision is kept fixed if EInit == EFnal. A NDO solver
	    which is *not* capable of handling approximate computation of the function
	    can ignore the values of these parameters.*/


			dblLastAlgParSG




		}; // end( dbl_par_type_SG )



	/*@} -----------------------------------------------------------------------*/
	/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/


	/** @name Friend classes
		@{ */

	/** The classes Stepsize and Deflection are "friends" of SubGradSolver. This is
	    done because Stepsize and Deflection objects may need some protected
	    data to work. An issue, however, is that derived classes from friend
	    classes are not friend, and therefore actual implementations of
	    Stepsize and Deflection cannot access data of SubGradSolver unless this
	    capability is explicitly provided by the base classes, who are friends.
	    This is why Stepsize and Deflection define a few methods that allow to
	    read protected data of SubGradSlver: so that any derived class can use them to
	    access to these data. These methods are, in fact, implemented at the end
	    of SubGradSolverc.Cpp. */

	  friend class Stepsize;
	  friend class Deflection;


	/*@} -----------------------------------------------------------------------*/




	/*--------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------*/
	/*---------------------------- CONSTRUCTOR ---------------------------------*/
	/*--------------------------------------------------------------------------*/
	/** @name Constructor
	    @{ */


	  /// constructor: does nothing special
	SubGradSolver();


	/*@} -----------------------------------------------------------------------*/






	/// destructor: it has to release all the Modifications
	virtual ~SubGradSolver();


	virtual ComputeConfig * get_ComputeConfig( bool all = false ,
					      ComputeConfig * ocfg = nullptr ) const override;

	virtual void set_ComputeConfig( ComputeConfig *scfg = nullptr ) override;











	/*@} -----------------------------------------------------------------------*/
	/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
	/*--------------------------------------------------------------------------*/
	/** @name Other initializations
	 *  @{ */


	/// set the (pointer to the) Block that the Solver has to solve

	 void set_Block( Block * block ) override;

	 /*--------------------------------------------------------------------------*/

	  virtual void set_par( const idx_type par , const int value ) override;

	 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	  virtual void set_par( const idx_type par , const double value ) override;


	 /*--------------------------------------------------------------------------*/



	/*--------------------------------------------------------------------------*/
	 /// set the ostream for the Solver log
	 /** Set the (pointer to) the ostream upon which the Solver is supposed to
	  * log all its algorithmic developments. Passing nullptr instructs the
	  * Solver to stay mum. The "verbosity" of such log is specified by
	  * intLogVerb [see set_par( int )], and setting intLogVerb = 0 should be
	  * equivalent to setting f_log = nullptr.
	  *
	  * The method is virtual because derived classes may need to do something to
	  * react to changes of the output log than storing the new log stream. */

	 void set_log( std::ostream *log_stream = nullptr ) override;








	 /*--------------------------------------------------------------------------*/

	/*@} -----------------------------------------------------------------------*/


	   /*--------------------------------------------------------------------------*/
	   /*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
	   /*--------------------------------------------------------------------------*/
	   /** @name Solving
	    *  @{ */


	   virtual int compute( bool changedvars = true ) override;




	/*@} -----------------------------------------------------------------------*/
	/*---------------------- METHODS FOR READING RESULTS -----------------------*/
	/*--------------------------------------------------------------------------*/
	/** @name Accessing the found solutions (if any)
	 *  @{ */


	 // the lb for the SubGradSolver is -Inf<double()> ???
	 virtual OFValue get_lb( void ) override;

	/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	 virtual OFValue get_ub( void )  override ;

	/*--------------------------------------------------------------------------*/

	 virtual bool has_var_solution( void ) override;
	  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	  virtual bool has_dual_solution( void ) override;

	 /*--------------------------------------------------------------------------*/
	 /*
	  virtual bool is_var_feasible( void ) override { return( true ); }

	  virtual bool is_dual_feasible( void ) override { return( true ); }
	 */
	 /*--------------------------------------------------------------------------*/
	  /// write the "current" solution
	  virtual void get_var_solution( Configuration *solc = nullptr ) override;

	 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	  /// write the "current" dual solution
	  virtual void get_dual_solution( Configuration *solc = nullptr ) override;

	 /*--------------------------------------------------------------------------*/

	  virtual bool new_var_solution( void ) override;

	 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	  virtual bool new_dual_solution( void )  override;
	 /*--------------------------------------------------------------------------*/
	 /*
	  virtual void set_unbounded_threshold( const OFValue thr ) override { }
	 */

	 /*--------------------------------------------------------------------------*/
//
//	  virtual bool has_var_direction( void ) override { return( true ); }
//
//	 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
//
//	  virtual bool has_dual_direction( void ) override { return( true ); }
//
//	 /*--------------------------------------------------------------------------*/
//	  /// write the current direction
//	  virtual void get_var_direction( Configuration *dirc = nullptr ) override
//	  {
//	   }
//
//	 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
//	  /// write the current dual direction
//	  virtual void get_dual_direction( Configuration *dirc = nullptr ) override
//	  {
//	   }

	 /*--------------------------------------------------------------------------*/
	 /*
	  virtual bool new_var_direction( void ) override { return( false ); }

	  virtual bool new_dual_direction( void ) override{ return( false ); }
	 */


/*@}*/








      /*--------------------------------------------------------------------------*/
	 /*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
	 /*--------------------------------------------------------------------------*/
	 /** @name Reading the data of the problem
	     @{ */





	 /*@}------------------------------------------------------------------------*/


	    /*@} -----------------------------------------------------------------------*/
	    /*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
	    /*--------------------------------------------------------------------------*/
	    /** @name Handling the parameters of the Solver
	     *  @{ */

	     virtual idx_type get_num_int_par( void ) const override
	     {
	      return( idx_type( intLastAlgParSG ) );
	      }

	    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	     virtual idx_type get_num_dbl_par( void ) const override
	     {
	      return( idx_type( dblLastAlgParSG ) );
	      }

	    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	     virtual idx_type get_num_str_par( void ) const override
	     {
	      return( idx_type( strLastAlgPar ) );
	      }

	     /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	     virtual idx_type dbl_par_str2idx( const std::string & name ) const override
	      {
	    	 // these may be many enough as to warrant using a map
	    	 const auto it = dbl_pars_map.find( name );
	    	 if( it != dbl_pars_map.end() )
	    	  return( it->second );
	    	 else
			   return( Solver::dbl_par_str2idx( name ) );
	}

	     /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
	      virtual idx_type int_par_str2idx( const std::string & name ) const override
	      {
	    	  // these may be many enough as to warrant using a map
	    	  const auto it = int_pars_map.find( name );
	    	  if( it != int_pars_map.end() )
	    	   return( it->second );
	    	  else

	       return( Solver::int_par_str2idx( name ) );
	       }

	     /*--------------------------------------------------------------------------*/

	      virtual const std::string & dbl_par_idx2str( const idx_type idx )
	        const override
	       {

	    	  	        if( ( idx >= dblSGPar2 ) && ( idx < dblLastAlgParSG ) )
	    	  	         return( dbl_pars_str[ idx - dblSGPar2 ] );
	    	  	        else
	    	  	         return( Solver::dbl_par_idx2str( idx ) );
	        }

	      /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/


	      /*--------------------------------------------------------------------------*/

	       virtual int get_int_par( const idx_type par ) const override;

	      /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	       virtual double get_dbl_par( const idx_type par ) const override;

	      /*--------------------------------------------------------------------------*/

	      /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	       virtual const std::string & int_par_idx2str( const idx_type idx )
	        const override
	       {

	        if( ( idx >= intSGPar1 ) && ( idx < intLastAlgParSG ) )
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
	        	if( ( par >= intSGPar1 ) && ( par < intLastAlgParSG ) )
	        	  return( dflt_dbl_par[ par - intSGPar1 ] );
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
	        	if( ( par >= dblSGPar2 ) && ( par < dblLastAlgParSG ) )
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

	        virtual const std::string & get_dflt_str_par( const idx_type par ) const override {
	         static const std::string empty = "";
	         return( empty );
	         }
	        /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	        /** get the configuration params for deflection */
	        Configuration* getConfdeflection() const;


	        /** set the configuration for deflection, the params configuration is a
	         * SimopleConfiguraztion<double> */
	        void setConfdeflection( Configuration*& confdeflection);


	        /** get the configuration for stepsize */
	        Configuration* getConfstepsize() const ;


	        /** set the configuration for stepsize, the configuration is a
	        	         * SimopleConfiguraztion<double> */
	        void setConfstepsize( Configuration*& confstepsize);

	      /*@} -----------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE TYPES -------------------------------*/
/*--------------------------------------------------------------------------*/


	 SMSpp_insert_in_factory_h;



/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

	 /*--------------------------------------------------------------------------*/

	 void Log1( void );
	 void Log2( void );
	 void Log3( void );
	 void Log4( void );

	/*--------------------------------------------------------------------------*/


	   void AllocateUBLB( const Variable::Index strt , const Variable::Index num );

	/*--------------------------------------------------------------------------*/

	   double lbnd( const Variable::Index i );

	/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

	  double  ubnd( const Variable::Index i );

	/*--------------------------------------------------------------------------*/

	  #if SubGrad_HANDLES_CONSTRAINTS
	   void EvalConst( void );
	  #endif

	   void EvalGi( Variable::Index wFi = Inf<Variable::Index>() );

	/*--------------------------------------------------------------------------*/

	   bool ProjFsb( std::vector<double> Lpt , bool & EmptySet );

	   void ProjTng( std::vector<double> Gpt , std::vector<double> Lpt ,
			 const Variable::Index strt = 0 , Variable::Index stp = Inf<Variable::Index>() );

	/*--------------------------------------------------------------------------*/

	   void InitStepsize( void );

	   void InitDeflection( void );

	/*--------------------------------------------------------------------------*/

	   void ChgLambdaHat( void );

	/*--------------------------------------------------------------------------*/

	   void UpdateSigma( void );

	/*--------------------------------------------------------------------------*/

	//   void UpdtLowerBound( void );

	/*--------------------------------------------------------------------------*/

	   //void MemDealloc( void );

	/*--------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------*/

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

	         void MemDealloc( void );

	      /*--------------------------------------------------------------------------*/



protected:

	   /*--------------------------------------------------------------------------*/
	   /*-------------------------- PROTECTED METHODS -----------------------------*/
	   /*--------------------------------------------------------------------------*/


	      void FiAndGi( const Variable::Index wFi = Inf<Variable::Index>() );

	     /**< Evaluates the function at the new point \f$\lambda_{i+1}\f$, i.e.,
	      * \f$f(\lambda_{i+1})\f$, and it either computes a subgradient
	      * \f$g_{i+1} \in \partial f(\lambda_{i+1})\f$, or, if the point is
	      * infeasible, a constraint. */

	     /*--------------------------------------------------------------------------*/

	         void GotoLambda1( void );

	      /**< Move the current point to \f$\lambda_{i+1}\f$. */

	      /*--------------------------------------------------------------------------*/
	   /*--------------------------------------------------------------------------*/

	   	 /**< The method is used within Solve(), and its job is to compute the
	   	  * direction \f$d_i\f$ that appears in the formula of \f$\lambda_{i+1}\f$.
	   	  * The direction is actually saved after the variables generation because i)
	   	  * the subgradient keeps in memory just the current direction, and ii) by
	   	  * generating/removing variables the direction may quickly come to be
	   	  * deteriorated. When the variables generation is ended [see GetFiStatus() in
	   	  * FiOracle.h], SaveDir() must be called in order to update the direction.
	   	  *
	   	  * In addition, the deflection coefficient is computed inside FormD(). As for
	   	  * the stepsize, the computation is performed within FormD() only if the
	   	  * scheme is <em>deflection-restricted</em>. */

	   	 void FormD( void ) ;

	   	 /*--------------------------------------------------------------------------*/

	   	 void SaveDir( void );

	   	 /**< The method is combined with FormD() [see above] to carry out the
	   	  * direction computation. Moreover, the stepsize is determined into SaveDir()
	   	  * when the adopted scheme is <em>stepsize-restricted</em>. */

	   	 /*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  ------------------*/
/*--------------------------------------------------------------------------*/


	   std::vector< C05Function * > v_c05f; /* the vector of the components of the
	                                           sum function */
	   LinearFunction * lf; ///< the 0-th component of the sum function



	    std::vector< std::pair <C05Function::LinearizationName , C05Function::LinearCombination > > zA;
	    /* the vector of the pairs  important linearization name and the
	       linear combination used to form it */

	   Vec_p_Var LamVcblr;    // the set of indices of Lambda

	   Variable::Index KSGPar1;       ///< projection-strategy parameters
//
	   double KSGPar2;       ///< incremental factor
//
	   Variable::Index KSGPar3;       ///< scheme: stepsize(deflection)-restricted
//
	   bool KSGPar4;        ///< control if \f$ \hat{\lambda}_i \f$ is kept
//
	   Variable::Index KSGPar5;       ///< seed
//
	    Configuration *confstepsize;
	    Configuration *confdeflection;



	   Deflection *deflection;   ///< pointer to the Deflection class

	   Stepsize *stepsize;       ///< pointer to the Stepsize class
//


//
//	   Variable::Index MaxNumVar;     ///< maximum number of variables
//
	   Function::Vec_FunctionValue LambdaBar;     ///< the stability center \f$\bar{\lambda}_i\f$
	   double FiBar;         ///< full function value at \a LambdaBar
//
	   Function::Vec_FunctionValue Lambda;        /**< the current point \f$ \lambda_i \f$
//	                             or the trial point \f$ \lambda_{i+1} \f$ */
	   double FiLambda;      ///< full function value at \a Lambda
//
	   Function::Vec_FunctionValue LambdaHat;     ///< the point \f$ \hat{\lambda}_i \f$
	   double FiHat;         ///< full function value at \a HLmb
//
//
	   bool LHasChgd;       /**< true if Lambda has changed since the latest call
//	                         * to FiAndGi() */
	   bool LHasProj;       /**< true if \a Lambda  has projected in the current
//	                         * iteration */
//
	   bool KpBstL;         ///< if LambdaBest has to be kept
	   double FiBest;        ///< the best value of \f$f\f$ found so far
	   Function::Vec_FunctionValue LambdaBest;    ///< the best point found so far
//
	   double LowerBound;    ///< Lower Bound over the full function \f$f\f$
	   bool TrueLB;         /**< true if LowerBound is a "true" lower bound rather
//	                         * than just the "minus infinity" */
//
	   Function::Vec_FunctionValue Gi;            ///< Gi[ wFi ]( Lambda ), the subgradient
	   Function::Vec_FunctionValue dir;           ///< the direction \f$ d_i \f$
//
	   double alpha;         ///< the deflection parameter \f$ \alpha_i \f$
	   double step;          ///< the stepsize \f$ \nu_i \f$
//
	   double Sigma;         /**< the linearization error \f$ \sigma_i\f$ of
//	                         * \f$ g_i\f$ at \f$\bar{\lambda}_{i}\f$ */
	   double Epsilon;       /**< the linearization error \f$\epsilon_i\f$ of
//	                         * \f$d_i\f$ at \f$\bar{\lambda}_i\f$ */
//
	   double SigmaHat;      /**< the linearization error \f$ \hat{\alpha}_i \f$ of
//	                         * \f$g_i\f$ with respect to \f$\hat{\lambda}_i\f$ */
	   double HatEpsilon;    /**< the linearization error \f$ \hat{\epsilon}_i \f$
//				 * of \f$d_i\f$ with respect to \f$\hat{\lambda}_i\f$
//				 */
//
	   std::vector<Variable::Index> SGBase;    ///< the set of indices of Gi[ wFi ]( Lambda )
//
	   double NrmGi;         ///< the (squared) subgradient's norm
	   double NrmDir;        ///< the (squared) direction's norm
//
	   double dGk;           ///< scalar product \f$ < d_i , g_i > \f$
	   double dM1Gk;         ///< scalar product \f$ < d_{i-1} , g_i > \f$
//
	   bool dM1GkDone;      /**< true if the scalar product
//				   \f$d_{i-1}^{\top} g_i\f$ has been computed */
//
	   Function::Vec_FunctionValue ub;            ///< upper bounds on the variables
	   Function::Vec_FunctionValue lb;            ///< lower bounds on the variables
//

//
	   bool ZeroComp;       ///< true if Fi() comes with the 0-th component
	   Variable::Index NItIncr;       /**< number of incremental iterations after an outer
//				   iteration */
	   bool InnIter;        ///< if true, the current iteration is an inner one
//
	   std::vector<Variable::Index> Seq;   /**< vector containing the randomly shuffled
//	                           components of the function \f$f\f$ */
	   bool DirPos;         /**< indicates where the direction \f$d_i\f$ in
//				   located in  the oracle */
	   Variable::Index MultBse;
//
	   Variable::Index CSSCntr;       ///< counter of consecutive SS
	   Variable::Index CNSCntr;       ///< counter of consecutive NS
//
	   Variable::Index CSmallStep;    ///< counter of consecutive \a short step
	   bool DoSS;           ///< SS vs NS
//
	   FiStatus fs;  ///< FiOracle status
	   //C05Function::




	  // C05Function* Oracle;  ///< (pointer to) the C05Function for Fi
//
	   bool EmptySet;       ///< true, if the feasible set is empty



	   Variable::Index NumVar;      ///< (current) number of variables

	   double EpsLin;      ///< optimality related parameter: relative precision

	   double tStar;       ///< optimality related parameter: "scaling" of Fi

	   Variable::Index ParIter;     ///< nuber of iterations in this run


	   Variable::Index MaxIter;     ///< maximum number of iterations
		double MaxTime;     ///< maximum time (in seconds) for each call to Solve()

		double EInit;       ///< precision-related parameter: initial precision
		double EFnal;       ///< precision-related parameter: final precision
		double EDcrs;       ///< precision-related parameter: rate of decrease
		int EStps;         ///< precision-related parameter: number of steps

//		NDOStatus Result;  ///< result of the latest call to Solve()


		//Variable::Index NrFi;        ///< number of components of Fi()

		Variable::Index SCalls;      ///< nuber of calls to Solve() (the current included)

		Variable::Index FiEvaltns;   ///< total number of Fi() calls
		Variable::Index GiEvaltns;   ///< total number of Gi() calls


		//std::ostream *NDOLog;   ///< the output stream object for log purposes
		//char NDOLLvl;      ///< the "level of verbosity" of the log

//		OPTtimers *NDOt;   ///< OPTtimer for timing purposes
		int Result;        ///< result of the latest call to Solve()

		double RelAcc; ///< relative accuracy for declaring a solution optimal
		double AbsAcc; ///< absolute accuracy for declaring a solution optimal
		double UpCutOff;    ///< upper cutoff for stopping the algorithm
		double LwCutOff;    ///< lower cutoff for stopping the algorithm
		double RAccSol;     ///< maximum relative error in any reported solution
		double AAccSol;     ///< maximum absolute error in any reported solution
		double FAccSol;     ///< maximum constraint violation in any reported solution

		int MaxSol; ///< maximum number of different solutions to report
		int LogVerb; ///< "verbosity" of the log



		/*--------------------------------------------------------------------------*/

		 const static std::vector<int> dflt_int_par;
		 ///< the (static const) vector of int parameters default values

		 const static std::vector<double> dflt_dbl_par;
		 ///< the (static const) vector of double parameters default values

		 const static std::vector< std::string > int_pars_str;
		 ///< the (static const) vector of int parameters names

		 const static std::vector< std::string > dbl_pars_str;
		 ///< the (static const) vector of double parameters names

		 const static std::map< std::string , idx_type > int_pars_map;
		  ///< the (static const) map for int parameters names

		 const static std::map< std::string , idx_type > dbl_pars_map;
		 ///< the (static const) map for double parameters names

};

#endif /* SUBGRADSOLVER_H_ */
