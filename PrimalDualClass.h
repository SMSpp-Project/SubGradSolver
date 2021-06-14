/*--------------------------------------------------------------------------*/
/*--------------------------- File PrimalDualClass.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * An intermediate class between Stepsize and Deflection and PrimalDual.
 * 
 * \version 0.1
 *
 * \date 21 - 05 - 2021
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
 * Copyright &copy 2001 - 2020 by Antonio Frangioni, Enrico Gorgone, Laura Fiorella.
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __PrimalDualClass
 #define __PrimalDualClass/* self-identification: #endif at the end of the file */


/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/
namespace SMSpp_di_unipi_it
{

 class SubGradSolver;     // forward declaration of class SubGradSolver

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Deflection ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup PrimalDualClass_CLASSES Classes in PrimalDualClass.h
    @{ */

/**
 * A class created to make factory system more efficient.
 */

class PrimalDualClass : public Deflection , public Stepsize
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
 
 using VarValue = Function::FunctionValue;
 using c_VarValue = Function::c_FunctionValue;

 using Vec_VarValue = Function::Vec_FunctionValue;
using Index = Function::Index;
using idx_type=ThinComputeInterface::idx_type;
/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{ */

/** Constructor of the class. It has no parameters apart the SubGradSolver object.*/

   PrimalDualClass( SubGradSolver * slvr=nullptr ): Deflection( slvr ) , Stepsize( slvr )
   {
    if( ! slvr )
    throw (std::invalid_argument( "PrimalDualClass: no subgradient solver " ));
    }  // end( PrimalDualClass::PrimalDualClass )

/**@}--------------------------------------------------------------------------*/

 static PrimalDualClass * new_PrimalDualClass( const std::string & classname, SubGradSolver * slvr = nullptr ) {
  const std::string classname_( SMSpp_classname_normalise(
					        std::string( classname ) ) );
  const auto it = PrimalDualClass::f_factory().find( classname_ );
  if( it == PrimalDualClass::f_factory().end() )
   throw( std::invalid_argument( classname +
				 " not present in PrimalDualClass factory" ) );
  return( ( it->second )( slvr ) );
  }
  
   const std::string & classname( void ) { return( private_name() ); }
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~PrimalDualClass() {}

/**@} ----------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:
 
/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/
/**
 * Factory definition
 **/
typedef boost::function< PrimalDualClass * (SubGradSolver *) > PrimalDualClassFactory;
// type of the factory of PrimalDualClass

typedef std::map< std::string, PrimalDualClassFactory > PrimalDualClassFactoryMap;
// type of the map between strings and the factory of PrimalDualClass

static PrimalDualClassFactoryMap & f_factory( void );

static void static_initialization(void){};
 
/*--------------------------------------------------------------------------*/
private:

virtual const std::string & private_name( void ) const = 0;

 };  // end( class PrimalDualClass )

/** @} end( group( PrimalDual_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

};  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif 

/*--------------------------------------------------------------------------*/
/*------------------------ End File PrimalDualClass.h ---------------------------*/
/*--------------------------------------------------------------------------*/
