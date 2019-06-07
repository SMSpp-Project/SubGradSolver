/*--------------------------------------------------------------------------*/
/*--------------------------- File Deflection.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the abstract class <tt>Deflection</tt>, which defines the
 * interface for computing the delflection coefficient within the SubGradSolver
 * class [see SubGradSolver.h].
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

#ifndef DEFLECTION_H_
 #define DEFLECTION_H_ /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS  ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_VOL 0

/* If LOG_VOL > 0, the Volume class produces a log of its activities on the
   ostream object and at the "level of verbosity" set with the method
   SetVOLLog() [see below]. */

#if( LOG_VOL )
 #define VOLLOG( l , x ) if( VOLLLvl > l ) *VOLLog << x
 #define VOLLOG2( l , c , x ) if( ( VOLLLvl > l ) && c ) *VOLLog << x
#else
 #define VOLLOG( l , x )
 #define VOLLOG2( l , c , x )
#endif




#include "SubGradSolver.h"
#include "Stepsize.h"
#include <algorithm>


/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


 class Stepsize;

 class SubGradSolver;     // forward declaration of class SubGradSolver



/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Deflection ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Deflection_CLASSES Classes in Deflection.h
    @{ */

/** The class Deflection provides an interface for the deflection rule (DR),
 *  to be used in the SubGradSolver solver [see SubGradSolver.h].
 *
 *  The aim of this class is to compute the deflection coefficient
 *  \f$ \alpha_i \f$, which is employed in the search direction formula as
 *  follows:
 *  \f[
 *   d_i = \alpha_i g_i + (1-\alpha_i)d_{i-1}\;.
 *  \f]
 *
 *  The user must extend the class to one or more deflection rules following
 *  this interface. */





class Deflection
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

/** Constructor of the class. It has no parameters apart the SubGradSolver object
 * using it because there is no "standard" deflection rule (except "no
 * deflection", i.e., \f$ \alpha_i = 1 \f$, which however in SubGradSolver is
 * treated by just setting the null Deflection object). */

   Deflection( SubGradSolver *slvr );

   static Deflection* new_Deflection( const std::string &classname , SubGradSolver *slvr  = nullptr );

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
   @{ */

/** The class outputs "log" information onto the ostream pointed by outs.
 * lvl controls the "level of verbosity" of the code: lvl == 0 means that
 * nothing at all is printed, and values larger than 0 mean increasing
 * amounts of information, the specific effect of each value being derived-
 * class-dependent. outs == 0 implies lvl == 0. */

   virtual void SetVOLLog( std::ostream *outs = 0 , const char lvl = 0 );

/*--------------------------------------------------------------------------*/
/** The method initializes the DR. It does nothing for the base class. */

   virtual void Format( void ) ;

/*@} -----------------------------------------------------------------------*/
/*-------------------- METHODS FOR DEFLECTION COMPUTATION ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
    @{ */

/** This method must be called before GetDFLCoeff() [see below].
 * This is indeed the core of every derived class, producing a new
 * deflection coefficient.
 *
 * Typically, the previous coefficient will be unavailable after the call to
 * NewDEF() [GetDFLCoeff()]. */

   virtual void NewDEF( void )=0;

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
    @{ */

/** Returns the deflection coefficient. This function must be called after
   NewDEF() [see above]. */

   virtual double GetDFLCoeff( void )=0;

/*--------------------------------------------------------------------------*/
/** Returns true if a <em>Serious Step</em> (SS) takes place. Otherwise a
 * <em>Null Step (NS)</em> occurs. By default, true is returned.
 *
 * Typically, a SS occurs when a \" good \" improvement of the function
 * \f$ f \f$ is obtained, usually it should be as good as the expected 
 * improvement [see Delta()]. */

   virtual const bool DoSS( void );

/*--------------------------------------------------------------------------*/
/** Returns the expected improvement in the objective function [see DoSS()].
 * This value is not useful for all DR. By default Inf<HpNum>() is returned.
 */

   virtual double Delta( void );

/*@}------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~Deflection() {}

/*@} -----------------------------------------------------------------------*/
   /*---------------------- PRIVATE PART OF THE CLASS -----------------------*/
   /*--------------------------------------------------------------------------*/

 private:



    virtual const std::string & private_name( void )  const = 0 ;



/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/





 protected:


    enum param_deflection{

    	last_param_deflection

    };


   /*--------------------------------------------------------------------------*/
      /*--------------------------- PROTECTED TYPES ------------------------------*/
      /*--------------------------------------------------------------------------*/

      typedef boost::function<Deflection*(SubGradSolver*)> DeflectionFactory;
       ///< type of the factory of Deflection

       typedef std::map<std::string,DeflectionFactory> DeflectionFactoryMap;
       ///< type of the map between strings and the factory of Deflection





/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods to read from SubGradSolver
 *
 * These methods are used to provide access to the information stored into
 * the SubGradSolver object. Deflection is a "friend" of SubGradSolver and therefore it
 * can read its protected data structures, but classes derived from
 * Deflection are not friend of SubGradSolver and cannot. This is why these methods
 * are defined in the base class (and implemented in SubGradSolver.C).
 * @{ */


       /// method incapsulating the Deflection factory
             /** This method returns the Deflection factory, which is a static object.
              * The rationale for using a method is that this is the "Construct On First
              * Use Idiom" that solves the "static initialization order problem". */

           static DeflectionFactoryMap & f_factory( void );


/*--------------------------------------------------------------------------*/
/** Returns the pointer to FiOracle. Thus, Deflection can ask directly the
 * FiOracle object for the function information. */

   //FiOracle* GetOracle( void );

/*--------------------------------------------------------------------------*/
/** Returns the stepsize \f$\nu_i\f$. */

   double GetStepsize( void );

/*--------------------------------------------------------------------------*/
/** Returns the norm of the SubGradSolverient \f$ g_i \f$. */

  const double GetGiNorm( void );

/*--------------------------------------------------------------------------*/
/** Returns the norm of the direction \f$ d_i \f$. */

  const double GetDNorm( void );

/*--------------------------------------------------------------------------*/
/** Returns the scalar product \f$ g_i^{\top} d_i \f$. */

  const double GetdGk( void );

/*--------------------------------------------------------------------------*/
/** Returns the linearization error \f$ \sigma_i \f$. */

  const double GetSigma( void );

/*--------------------------------------------------------------------------*/
/** Returns the linearization error \f$ \epsilon_i \f$. */

  const  double GetEpsilon( void );

/*--------------------------------------------------------------------------*/
/** Returns the full function \f$ f(\lambda_i) \f$. */

  const  double ReadFVal( void );


  /*--------------------------------------------------------------------------*/
  /** Returns the pointer to FiOracle. Thus, Deflection can ask directly the
   * FiOracle object for the function information. */

     std::vector<C05Function*> GetOracle( void );

  /*--------------------------------------------------------------------------*/

/*@} -----------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

  SubGradSolver *Solver;    ///< (pointer to) the SubGradSolver solver
  std::ostream *VOLLog;    ///< the output stream object
  char VOLLLvl;       ///< the "level of verbosity"

  SimpleConfiguration<std::vector< double > > *deflectionconf;

/*--------------------------------------------------------------------------*/

 };  // end( class Deflection )

/* @} end( group( Deflection_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/





class PrimalDual :  public Deflection, public Stepsize
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/


 public:



	using Deflection::f_factory;
	//using Stepsize::f_factory;

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

/** Constructor of the class, which derives from both Deflection and
   Stepsize. Since the constructor of PrimalDual is executed <em>after</em>
   the one of Stepsize, the following parameters specific for PrimalDual have
   to be found in the stream \a after those of the base class Stepsize [see the
   comments to the constructor of Stepsize]:

     -# bool average [true] if true simple averages are used, otherwise
                            weighted averages are used

   For this class, the general parameter LpsFct -defined in Stepsize.h-
   represents the gamma factor \a F used to set \f$\omega_i\f$. */



	enum param_deflectionPrimalDual{

		par_average = last_param_deflection , ///[true] if true simple averages are used, otherwiseweighted averages are used



	   	 last_param_deflectionPrimalDual

	    };


PrimalDual( SubGradSolver *slvr  );



/*@} -----------------------------------------------------------------------*/
 /*-------------------------- OTHER INITIALIZATIONS -------------------------*/
 /*--------------------------------------------------------------------------*/
 /** @name Other initializations
   @{ */

   void SetVOLLog( std::ostream *outs = 0 , const char lvl = 0 ) override;

   void Format( void ) override;

/*@} -----------------------------------------------------------------------*/
/*-----------------------  METHODS FOR DR and SR COMPUTATION  --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
     @{ */

/** The variant requires thats DR must be computed before SR. Hence, the call
    to NewDEF comes before the one to NewStep. */

   void NewDEF( void ) override;

   void NewStep( void ) override;

/*--------------------------------------------------------------------------*/

   double GetDFLCoeff( void ) override;

   double GetStepsize( bool StepIsIncr = false ) override;

/*--------------------------------------------------------------------------*/

  double GetLev( void ) override;


/*--------------------------------------------------------------------------*/

  double GetBeta( void ) override;

/*--------------------------------------------------------------------------*/

  void SetMaxBeta( const double alpha ) override;

/*@} -----------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

  SMSpp_insert_in_factory_h;

   void UpdateGamma( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   double delta;
   double Tmpdelta;
   double vk;

   double StepSize;
   double gamma;

   double beta1;
   bool average;

   double LipCnst;      // Lipschitz constant

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class PrimalDual )



/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Definition of the class Volume. This class implements the Volume algorithm.
    The method is revisited and incorporated in the SubGradSolver solver. The
    deflection coefficient is found solving a quadratic problem that
    involves the subgradient \f$ g_i \f$, the direction \f$ d_{i-1} \f$, and
    their respective linearization errors \f$ \sigma_i \f$  and
    \f$ \epsilon_{i-1} \f$ at the stability center \f$ \bar{\lambda}_i\f$,
    namely
    \f[
     \tau_i =
     \arg\min \left\{ \; \nu_{i-1} \left\| \tau g_i
                         + (1 - \tau) d_{i-1} \right\|^2/2
                         +  \tau \sigma_i(\bar{\lambda}_i)
                         + (1 - \tau)\epsilon_{i-1}(\bar{\lambda}_i)
        \;:\; \tau \in [0, 1]\; \right\}
    \f]  */

class Volume : public Deflection
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

/** Constructor of the class. The parameter `iStrm', if provided, is taken as
   a pointer to a istream from which the algorithmic parameters for the
   Volume are sequentially read in the following order. Each parameter must be
   placed at the beginning of a separate line, max 255 characters long, with
   all the rest of the line up to the first newline character '\n' (apart from
   a separating whitespace) being available for comments. Any line whose first
   character is '#' and any blank line is ignored. If 0 is passed, the file
   ends before reaching a given parameter, or some parameter is in the wrong
   format, each non-specified parameter is given a default value, shown in []
   below.

    -#   double tauInit   [0.05]  \f$ \tau_0 \f$

    -#   double tauMin    [1e-4]  safety threshold of \f$ \tau:~ \tau_{\min}\f$

    -#   double tauFactor [0.5]   factor \f$ \tau_{f} \f$

    -#   double tauIter   [50]    $\tau_p$

    -#   double m         [1e-3]  descent parameter */

  Volume( SubGradSolver *slvr );

 /*@} -----------------------------------------------------------------------*/
 /*-------------------------- OTHER INITIALIZATIONS -------------------------*/
 /*--------------------------------------------------------------------------*/
 /** @name Other initializations
   @{ */

   void SetVOLLog( std::ostream *outs = 0 , const char lvl = 0 ) override;

/*--------------------------------------------------------------------------*/

   void Format( void ) override;

/*@} -----------------------------------------------------------------------*/
/*-------------------- METHODS FOR DEFLECTION COMPUTATION ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
     @{ */

   void NewDEF( void ) override;

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
    @{ */

   double GetDFLCoeff( void ) override;

   const bool DoSS( void )override;

   double Delta( void ) override;

/*@} -----------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:


   enum param_deflectionVolume{

	   par_tauInit = last_param_deflection , ///[0.05]  \f$ \tau_0 \f$

	   par_tauMin, ///< [1e-4]  safety threshold of \f$ \tau:~ \tau_{\min}\f$

	   par_tauFactor, ///< [0.5]   factor \f$ \tau_{f} \f$

   	 par_tauIter, ///< [50] $\tau_p$

	 par_m, ///< [1e-3]  descent parameter

   	 last_param_deflectionVolume

    };


   SMSpp_insert_in_factory_h;
/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

   void Solve( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   double tauInit;     // initial tau
   double tauMin;      // the minimum value of tau
   double tauFactor;   // the decreasing factor
   Variable::Index tauIter;     // the number of iterations after a decrease
   double m;           // descent parameter

   double tau;         // the safeguarded parameter tau
   double lastvalue;   // last function value Fi()
   double alpha;       // deflection coefficient

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class Volume )



/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Deflection.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File Deflection.h ---------------------------*/
/*--------------------------------------------------------------------------*/
