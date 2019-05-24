/*--------------------------------------------------------------------------*/
/*--------------------------- File Stepsize.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Definition of the abstract <tt>Stepsize</tt> class, which defines the
 * interface for computing the stepsize rule (SR) used by the
 * <tt>SUBGRADSolver</tt> class [see SUBGRADSolver.h].
 *
 * \version 1.01
 *
 * \date 27 - 01 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Department of Informatics  \n
 *         University of Pisa \n \n
 *
 * \author Enrico Gorgone \n
 *         Department of Mathematics and Informatics \n
 *         University of Cagliari \n\n
 *
 * Copyright &copy 2001 - 2019 by Antonio Frangioni, Enrico Gorgone.
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef STEPSIZE_H_
 #define STEPSIZE_H_ /* self-identification: #endif at the end of the file */




#define LOG_STP 0

/* If LOG_STP > 0, the FumeroTV class produces a log of its activities on the
   ostream object and at the "level of verbosity" set with the method
   SetSTPLog() [see below]. */

#if( LOG_STP )
 #define STPLOG( l , x ) if( STPLLvl > l ) *STPLog << x
 #define STPLOG2( l , c , x ) if( ( STPLLvl > l ) && c ) *STPLog << x
#else
 #define STPLOG( l , x )
 #define STPLOG2( l , c , x )
#endif







/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SUBGRADSolver.h"
#include "C05Function.h"


#include <algorithm>





/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

 class SUBGRADSolver;     // forward declaration of class SUBGRADSolver


/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Stepsize --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Stepsize_CLASSES Classes in Stepsize.h
    @{ */

/** The class Stepsize provides an interface for the stepsize rule, to be
 *  used in the SUBGRADSolver solver [see SUBGRADSolver.h].
 *
 *  The main functions of the class are NewStep() and GetStepsize(). While
 *  NewStep() sets up a new step, the value \f$\nu_i\f$ is passed to the
 *  SubGrad solver by means of a call to GetStepsize().
 *
 *  The user must extend the class to one or more stepsize rules following
 *  this interface.  */

class Stepsize
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
	@{ */

/** Constructor of the class. The parameter `iStrm', if provided, is taken
 * as a pointer to a istream from which the algorithmic parameters for the
 * stepsize formula are sequentially read in the following order. Each must
 * be placed at the beginning of a separate line, max 255 characters long,
 * with all the rest of the line up to the first newline character '\n'
 * (apart from a separating whitespace) being available for comments. Any
 * line whose first character is '#' and any blank line is ignored.
 * If 0 is passed, the file ends before reaching a given parameter, or
 * some parameter is in the wrong format, each non-specified parameter is
 * given a default value, shown in [] below.
 *
 * Only a few parameters are deemed to be general enough to be put in the
 * basic Stepsize class; however, the same stream can be used for passing
 * the (possibly many more) algorithmic parameters to the derived classes.
 * Since the constructor of Stepsize is executed before the ones of the
 * derived classes, the specific parameters for the derived classes have just
 * to be found in the stream after those of the base class.
 *
 * HpNum LpsFct [1] a general scaling factor (for instance, it is the 
 *                  Lipschitz factor in the target value SR)  */




 Stepsize( SUBGRADSolver *slvr  ) ;// @suppress("Class members should be properly initialized")




 /*--------------------------------------------------------------------------*/
  /// construct a :Stepsize of specific type using the Stepsize factory
  /** Use the Stepsize factory to construct a :Stepsize object of type specified by
   * classname (a std::string with the name of the class inside). If there is
   * no class with the given name, exception is thrown. The optional parameter
   * father is the pointer to the father SUBGRADsolver, that (if provided) is passed
   * to the :Stepsize constructor inside the factory.
   *
   * Note that the method is static because the factory is static, hence it is
   * to be called as
   *
   *       Stepsize *Stepsize = Stepsize::newStepsize( someclass );
   *
   *
   * Note that the :Stepsize returned my this method is "empty": it contains no
   * (instance) data, and therefore it has to be explicitly initialized with
   * any of the corresponding methods (operator>>, serialize(), anything that
   * the specific :Stepsize class provides) before it can be used.
   *
   * For this to work, each :Stepsize has to:
   *
   * - add the line
   *
   *       SMSpp_insert_in_factory_h;
   *
   *   to its definition (typically, in the private part in its .h file);
   *
   * - add the line
   *
   *       SMSpp_insert_in_factory_cpp_1( name_of_the_class );
   *
   *   to exactly *one* .cpp file, typically that :Stepsize .cpp file. */



 static Stepsize* new_Stepsize( const std::string &classname , SUBGRADSolver *slvr  = nullptr );





/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
   @{ */

/** The class outputs "log" information onto the ostream pointed by outs.
   lvl controls the "level of verbosity" of the code: lvl == 0 means that
   nothing at all is printed, and values larger than 0 mean increasing
   amounts of information, the specific effect of each value being derived-
   class-dependent. outs == 0 implies lvl == 0. */

   virtual void SetSTPLog( std::ostream *outs = 0 , const char lvl = 0 );


/*--------------------------------------------------------------------------*/
/** The method initializes the stepsize rule. */

   virtual void Format( void );

/*@}------------------------------------------------------------------------*/
/*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Stepsize computation
    @{ */

/** This method must be called before GetStepsize() [see below]. This is
   indeed the core of every derived class, producing a new step.

   If the derived class implements a target value stepsize rule, the method
   must compute both the scalar \f$ \beta_i \f$ and the level
   \f$ f^{lev}_i \f$.

   Typically, the previous step will be unavailable after the call to
   NewStep() [see GetStepsize()].  */

   virtual void NewStep( void )=0;

/*--------------------------------------------------------------------------*/
/** Returns true if the stepsize requires the scalar product
   \f$ d_{i-1}^{\top} g_i\f$. The default implementation don't require this
   value. */

   virtual bool NeedsdkM1Gk( void ) ;


/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
    @{ */

/** This function must be called after NewStep() [see above ] and returns the
 * step value. StepIsIncr controls the stepsize formula in the incremental
 * version.
 *
 * The base class provides the step value for target stepsize rules, by
 * exploiting \f$ beta_k \f$ and \f$ f^{lev}_k \f$ computed by NewStep().
 * While for non-incremental step (StepIsIncr = false) for \f$ k = 1, 2 ,
 * \ldots \f$ the formula is
 * \f[
 *   \nu_k = \beta_k * ( f_k - f^{lev}_k ) / \| d_k \|^2\,,
 * \f]
 * for incremental step (StepIsIncr = true) for \f$ k = p(\Pi+1) , \ldots ,
 * (p+1)(\Pi+1) - 1 , p = 0, 1 , 2 \ldots \f$ the stepsize rule is
 * \f[
 *   \nu_k = \beta_k \frac{ f_{ p(\Pi + 1 ) } - f^{lev}_{ p( \Pi + 1 ) } }
 *                        { \chi( \Pi + 1 ) C^2 }
 * \f]
 * where \f$ \chi \f$ is the LpsFct factor and \f$ C \f$ is the Lipschitz
 * constant.
 *
 * Note that GetStepsize() only returns the value, the actual step
 * computation being done by NewStep(). */

   virtual double GetStepsize( bool StepIsIncr = false );

    
/*--------------------------------------------------------------------------*/
/** This function must be called after NewStep() [see above].
 *
 * The default implementation returns the level \f$ f^{lev}_i \f$ for a
 * target value stepsize rule [see below]. It should throw exception for
 * other forms of the stepsize rule. */

   virtual double GetLev( void );


/*--------------------------------------------------------------------------*/
/** This function must be called after NewStep() [see above].
 *
 * The default implementation returns the scalar \f$ \beta_i \f$ for a
 * target value stepsize rule [see below]. It should throw exception for
 * other forms of the stepsize rule. */

   virtual double GetBeta( void );

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
    @{ */

/** Changes the <em>maximum value</em> of \f$ \beta_i \f$ for a <em>target
 * value stepsize rule</em>.
 *
 * Typically, the method implements the safe rule of the stepsize-restricted
 * approach. Passing the deflection coefficient \f$ \alpha_i \f$ does exactly
 * this job. */

   virtual void SetMaxBeta( const double alpha );

/*@}------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~Stepsize() {};



 private:
   /*--------------------------------------------------------------------------*/
   /*-------------------------- PRIVATE METHODS -------------------------------*/
   /*--------------------------------------------------------------------------*/
    // Definition of Solver::private_name() (pure virtual)

    //virtual const std::string & private_name( void ) const = 0;

   virtual const std::string & private_name( void )  const = 0 ;
   //SMSpp_insert_in_factory_h;






/*@} -----------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

  protected:

   /*--------------------------------------------------------------------------*/
   /*--------------------------- PROTECTED TYPES ------------------------------*/
   /*--------------------------------------------------------------------------*/

   typedef boost::function<Stepsize*(SUBGRADSolver*)> StepsizeFactory;
    ///< type of the factory of StepSize

    typedef std::map<std::string,StepsizeFactory> StepsizeFactoryMap;
    ///< type of the map between strings and the factory of Stepsize

   /*--------------------------------------------------------------------------*/
   /*-------------------------- PROTECTED METHODS -----------------------------*/
   /*--------------------------------------------------------------------------*/

    /// method incapsulating the Stepsize factory
    /** This method returns the Stepsize factory, which is a static object.
     * The rationale for using a method is that this is the "Construct On First
     * Use Idiom" that solves the "static initialization order problem". */

    static StepsizeFactoryMap & f_factory( void );


/** It tries to provide a new target level \f$ f^{lev}_i \f$. It true is
 * returned, the level has been changed. By default, it returns false. */

   virtual bool UpdateTargetLevel( void );


/*--------------------------------------------------------------------------*/
/** @name Protected methods to read from Subgrad
 *
 * These methods are used to provide access to the information stored into
 * the SubGrad object. Stepsize is a "friend" of SubGrad and therefore it
 * can read its protected data structures, but classes derived from
 * Stepsize are not friend of SubGrad and cannot. This is why these methods
 * are defined in the base class (and implemented in SubGrad.C).
 * @{ */

/*--------------------------------------------------------------------------*/
/** Returns the deflection coefficient \f$\alpha_i\f$. */

   double GetCoeffDefl( void );

/*--------------------------------------------------------------------------*/
/** Returns the pointer to C05Function. Thus, stepsize can ask directly the
 * C05Function object for the function information. */

   C05Function* GetOracle( void );

/*--------------------------------------------------------------------------*/
/** Returns the norm of the subgradient \f$ g_i \f$. */

   const double GetGiNorm( void );

/*--------------------------------------------------------------------------*/
/** Returns the norm of the direction \f$ d_i \f$. */

   const double GetDNorm( void );

/*--------------------------------------------------------------------------*/
/** Returns the scalar product \f$ g_i^{\top} d_i\f$. */

   const double  GetdGk( void );

/*--------------------------------------------------------------------------*/
/** Returns the scalar product \f$ g_i^{\top} d_{i-1}\f$. */

   double  GetdkM1Gk( void );

/*--------------------------------------------------------------------------*/
/** Returns NItIncr, the parameter for incremental iterations [see SubGrad.h].
 */

   Variable::Index GetNItIcr( void );

/*--------------------------------------------------------------------------*/
/** Returns FiLambda [ see SubGrad.h ]. */

   const double ReadFkVal( void );
    
/*--------------------------------------------------------------------------*/
/** Returns FiBar [ see SubGrad.h ]. */
    
   const double ReadFiBar( void );

/*@}------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

enum param_stepsize{

	par_LpsFct, ///< scaling factor

	last_param_stepsize

};

  
  SUBGRADSolver *Solver;  ///< (pointer to) the SubGrad solver
  std::ostream *STPLog;  ///< the output stream object
  char STPLLvl;     ///< the "level of verbosity"

  double FiLev;      ///< the target level \f$ f^{lev}_i \f$
  double Beta;       ///< beta factor \f$ \beta_i \f$
  double MaxBeta;    ///< maximum value for beta factor

  double LpsFct;     ///< scaling factor


  SimpleConfiguration<std::vector< double > > *stepsizeconf;



 }; // end( class Stepsize )

/* @} end( group( Stepsize_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/






/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Definition of the class ColorTV. This class implements the target stepsize
    used in the original Volume algorithm. The method is based on classifying
    the iterations based on the obtained improvement of the function Fi(). */


class ColorTV : public Stepsize
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
    The ColorTV defines one main public types: */

 enum condition { green ,  ///< good step
		  yellow , ///< not bad, not good step
		  red      ///< bad step
                  };



 enum param_stepsizeColor{

 	 par_BetaZero = last_param_stepsize ,

	 par_greentestinvl, ///< maximum number of consecutive green iteration

	 par_yellowtestinvl, ///< maximum number of consecutive yellow iteration

	 par_redtestinvl, ///< maximum number of consecutive red iteration

	 last_param_stepsizeColor

 };


/*@}------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
	@{ */

/**  Constructor of the class. Since the constructor of ColorTV is executed
   after the one of Stepsize, the following parameters specific for ColorTV
   have to be found in the stream <em>after</em> those of the base class
   [see the comments to the constructor of Stepsize]:

     -# double BetaZero      [1] initial value of beta

     -# Index greentestinvl [1] how many consecutive green iterations
                                are allowed before increasing beta

     -# Index yellowtestinvl [400] how many consecutive yellow iterations
                             are allowed before decreasing beta

     -# Index redtestinvl    [10] how many consecutive red iterations are
                             allowed before decreasing beta */

  ColorTV( SUBGRADSolver *slvr );


/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
   @{ */

   void SetSTPLog( std::ostream *outs = 0 , const char lvl = 0 ) override;

/*--------------------------------------------------------------------------*/

   void Format( void ) override;

/*@}------------------------------------------------------------------------*/
/*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Computing the stepsize
     @{ */

    void NewStep( void ) override;

/*--------------------------------------------------------------------------*/

   bool NeedsdkM1Gk( void ) override { return( true );
    }

/*@} -----------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

  protected:

    bool UpdateTargetLevel( void ) override;

/*@}------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

    /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE METHODS -------------------------------*/
    /*--------------------------------------------------------------------------*/

     //SMSpp_insert_in_factory_h;

    /*---------------------------- PRIVATE FIELDS ------------------------------*/

     /* manual expansion of SMSpp_insert_in_factory_h to avoid including the
      * whole of SMSTypedefs.h. */
     static class _init {
     public:
      _init( void );
      } _initializer;

     virtual const std::string & private_name( void ) const override;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

    const condition coloring( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   Variable::Index greentestinvl;   ///< maximum number of consecutive green iteration
   Variable::Index yellowtestinvl;  ///< maximum number of consecutive yellow iteration
   Variable::Index redtestinvl;     ///< maximum number of consecutive red iteration

   Variable::Index lastgreeniter;   ///< the last green iteration
   Variable::Index lastyellowiter;  ///< the last yellow iteration
   Variable::Index lastrediter;     ///< the last red iteration

   double BetaZero;
   double LwrBnd;     ///< lower bound
   double lastvalue;  ///< the previous function value

   double FiLmbd;     ///< current FiLambda
   Variable::Index NrIter;

/*--------------------------------------------------------------------------*/

 };  // end( class ColorTV )






/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES CLASS FumeroTV------------------*/
/*--------------------------------------------------------------------------*/
/** Definition of the class FumeroTV. This class implements a target value
    stepsize. At the beginning \f$ f^{lev} \f$ is set to the given lower bound
    and exponentially goes towards \f$ f^{rec} \f$. The method works in two
    phases. In the first one \f$ \beta_i \f$ is decremented, while in the
    second one \f$ \beta_i \f$ is still decremented but also incremented
    after a pre-set number of consecutive improving iterations. */

class FumeroTV : public Stepsize
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

	enum param_stepsizeFumera{

	 	 par_sigmaMin = last_param_stepsize , ///< [ 1e-4 ] tolerance of the function \f$ sigma_r \f$

		 par_r1, ///< [ 2 ]    sigma_{r1} = 1/2

		 par_BetaZero, ///< [ 1 ]    initial value of beta

		 par_eta1, ///< [ 10 ]   threshold on number of failures (phase I)

		 par_eta2, ///< threshold on number of failures (phase II)

		 last_param_stepsizeFumera

	 };


/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{ */

/** Constructor of the class. Since the constructor of FumeroTV is executed
   after the one of Stepsize, the following parameters specific for the
   FumeroTV have to be found in the stream <em>after</em> those of the base
   class [see the comments to the constructor of Stepsize]:

     -# double sigmaMin  [ 1e-4 ] tolerance of the function \f$ sigma_r \f$

     -# Index r1        [ 2 ]    sigma_{r1} = 1/2

     -# double beta0     [ 1 ]    initial value of beta

     -# Index eta1      [ 10 ]   threshold on number of failures (phase I)

     -# index eta2      [ 10 ]   threshold on number of failures (phase II)
     */

 FumeroTV( SUBGRADSolver *slvr  );


/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
   @{ */

   void SetSTPLog( std::ostream *outs = 0 , const char lvl = 0 ) override;

/*--------------------------------------------------------------------------*/

   void Format( void ) override;

/*@}------------------------------------------------------------------------*/
/*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Computing the stepsize
    @{ */

    void NewStep( void ) override;

/*@} -----------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

     bool UpdateTargetLevel( void ) override;

/*@}------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/


    SMSpp_insert_in_factory_h;


    double ExpFun( Variable::Index ri );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   double FiLambda;   ///< FiLambda

   double LwrBnd;     ///< lower bound
   double FiRef;      ///< the best solution found so far

   double sigma;      ///< current value of \f$ \sigma(r)\ f$
   double sigmaMin;   ///< minimum threshold of \f$ \sigma(r) \f$

   double BetaZero;

   Variable::Index rc;         ///< current value of r
   Variable::Index r1;         ///<
   Variable::Index r2;

   Variable::Index etac;       ///< last improvement of FiLev
   Variable::Index eta1;       ///< thresholds for max number of failures
   Variable::Index eta2;

 }; // end( class FumeroTV )






/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES CLASS Poljak--------------------*/
/*--------------------------------------------------------------------------*/
/** Definition of the class Polyak. This class implements  a target value
    stepsize rule, whereby \f$\beta_i\f$ and \f$f^{lev}_i\f$ are constant
    (they do not depend on $i$).                                          */

class Polyak : public Stepsize
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

	enum param_stepsizePolyak{

		 	 par_Beta = last_param_stepsize , ///< [ 1e-2 ] beta value.

			 last_param_stepsizePolyak

		 };




/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
   @{ */

 Polyak( SUBGRADSolver *slvr ) ;

 /**<  Constructor of the class. Since the constructor of Polyak is executed
    after the one of Stepsize, the following parameters specific for the
    Polyak have to be found in the stream after those of the base class [see
    the comments in the constructor of Stepsize].

       -# HpNum Beta [ 1e-2 ] beta value. */

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
   @{ */

   void SetSTPLog( std::ostream *outs = 0 , const char lvl = 0 ) override;

/*--------------------------------------------------------------------------*/

   void Format( void ) override;

 /*@}------------------------------------------------------------------------*/
 /*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
 /*--------------------------------------------------------------------------*/
 /** @name Computing the stepsize
    @{ */

    void NewStep( void ) override;

/*@} -----------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

     bool UpdateTargetLevel( void ) override;

/*@}------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:



    SMSpp_insert_in_factory_h;
/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class Polyak )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/





#endif  /* Stepsize.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File Stepsize.h -----------------------------*/
/*--------------------------------------------------------------------------*/
