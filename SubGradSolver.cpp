/*--------------------------------------------------------------------------*/
/*------------------------ File SubGradSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
* Implementation of the SubGradSolver class.
*
* \version 0.1
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
* Copyright &copy 2019 - 2020 by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/
#include "SubGradSolver.h"

#include "Deflection.h"
#include "Stepsize.h"

#include "PrimalDual.h"
#include "FumeroTV.h"
#include "ColorTV.h"
#include "Polyak.h"
#include "Volume.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define VLOG( l , x ) if( LogVerb > l ) *f_log << x 
#define VLOG2( l , c , x ) if( ( LogVerb > l ) && c ) *f_log << x

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/
Index myrandom( Index i ) { return( rand() % i ); }

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SubGradSolver to the Solver factory
SMSpp_insert_in_factory_cpp_0( SubGradSolver ); 

SMSpp_insert_in_factory_cpp_1( (Volume ) ); 

SMSpp_insert_in_factory_cpp_1( (FumeroTV) ); 

SMSpp_insert_in_factory_cpp_1( ColorTV ); 

SMSpp_insert_in_factory_cpp_1( Polyak ); 

const std::string & PrimalDual::_private_name( void ) {
static const std::string _name( "PrimalDual" );
 return( _name );
 }

const std::string & PrimalDual::private_name( void ) const {
return( PrimalDual::_private_name() );
 }

PrimalDual::_init::_init( void ) {
 Stepsize::f_factory()[ PrimalDual::_private_name() ] =
                         boost::factory< PrimalDual * >();
 Deflection::f_factory()[ PrimalDual::_private_name() ] =
                         boost::factory< PrimalDual * >();
 Deflection::static_initialization();
 Stepsize::static_initialization();
 };

PrimalDual::_init PrimalDual::_initializer;


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

// define and initialize here the vector of int parameters names
const std::vector< std::string > SubGradSolver::int_pars_str =
			{ "intSGPar1" , "intSGPar3" , "intSGPar4","intSGPar5", "intRstAlgPar" };

// define and initialize here the vector of double parameters names
const std::vector< std::string > SubGradSolver::dbl_pars_str =
			{ "dblSGPar2" ,"dbltStar", "dblEpsLin","dblEInit","dblEFnalt","dblEDcrs","dblEStps" };			

// define and initialize here the map for int parameters names
const std::map< std::string , SubGradSolver::idx_type > SubGradSolver::int_pars_map =
                   { { "intSGPar1"  , SubGradSolver::intSGPar1  } ,
		     { "intSGPar3" , SubGradSolver::intSGPar3 } ,
		     { "intSGPar4" , SubGradSolver::intSGPar4 } ,
		     { "intSGPar5" , SubGradSolver::intSGPar5 },
			 {"intRstAlgPar", SubGradSolver::intRstAlgPar}};

// define and initialize here the map for double parameters names
const std::map< std::string , SubGradSolver::idx_type > SubGradSolver::dbl_pars_map =
                   { { "dblSGPar2" , SubGradSolver::dblSGPar2 } ,
		     { "dbltStar" , SubGradSolver::dbltStar } ,
			 { "dblEInit" , SubGradSolver::dblEInit } ,
			 { "dblEFnalt" , SubGradSolver::dblEFnalt } ,
			 { "dblEDcrs" , SubGradSolver::dblEDcrs } ,
			 { "dblEStps" , SubGradSolver::dblEStps } };


// define and initialize here the default int parameters
const std::vector<int> SubGradSolver::dflt_int_par =
        {    10 ,  // intSGPar1
			100 ,  // intSGPar3
			  0 ,  // intSGPar4
			  0 ,  // intSGPar5
			  2	   // intRstAlgPar
			 };

// define and initialize here the default double parameters
const std::vector<double> SubGradSolver::dflt_dbl_par =
           { 1,        //dblSGPar2
        	 1e2 ,    // dbltStar
			 1e-2 ,   // dblEpsLin
			 1e6 ,    // dblEInit
			 0.95 ,   // dblEFnalt
			 - 1 ,    // dblEDcrs
			 - 1     // dblEStps
               };

/*-------------------------------------------------------------------------*/

static const double Multplier = 1.0;

static const unsigned char PrjSG =  1;   // project of Gi[wFi] = \f$g_i\f$
static const unsigned char PrjDIRM1 = 2; // project of DirM1 = \f$d_{i-1}\f$
static const unsigned char PrjDIR = 4;   // project of Dir = \f$d_i\f$

static const unsigned char SafeRULE = 1;  // safe rule
static const unsigned char StepRST = 2;   // Stepsize-restricted

static const unsigned char RstAlg =  1;  // don't reset algorithmic parameters
static const unsigned char RstCrr =  2;  // don't reset current point
static const unsigned char RstSbg =  4;  // don't reset subgradients
//static const unsigned char RstCnt =  8;  // don't reset constraints
static const unsigned char RstFiV = 8;  // don't reset FiVals

//static const bool first = true;          // indicate the position of the
//static const bool second = false;        // direction

 /*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF SubGradSolver --------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
 
void SubGradSolver::set_Block( Block * block ){
 if( f_Block ) {  // changing from a previous oracle - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  guts_of_destructor();   // deallocate memory
  }
 Solver::set_Block( block );  // attach to the new Block

 if( ! f_Block )  // that was actually clearing the Block
  return;         // all done

 /* Two types of block can be handled by the BundleSolver:

     1. Only one single non-smooth function
     2. a sum of some non-smooth functions

    The algorithm here developed aims at solving non-constrained non-smooth
    optimization. The block can have box constraints at the most.
    It is expected the block to have in the first case a FRealObjective
    whose the function is a C05Function one and having no children, while in
    the second case a FRealObjective one whose the function is a
    LinearFunction and having as many sub-blocks as the number of components.
    In the latter case, each sub-block must not contain any Variable or
    Constraint.cVariable may have a lower and upper bound. If the lower
    bound  has a finite value, it must be 0. */

 const auto & sb = f_Block->get_nested_Blocks();

 if( sb.empty() ) {
  // the objective function of the block must be a C05Function  - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
  if( ! obj )
   throw( std::logic_error( "objective is not a FRealObjective" ) );

  auto c05f = dynamic_cast< C05Function * >( (obj)->get_function() );
  if( ! c05f )
   throw( std::logic_error( "the objective is not a C05Function" ) );

  v_c05f.push_back( c05f );
  f_lf = nullptr;
  }
 else {
  // the objective function of the block must be a LinearFunction- - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ! f_Block->get_objective() )  // there is no Objective
   f_lf =  nullptr;
  else {
   auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
   if( ! obj )
    throw( std::logic_error( "the objective is not a real function" ) );

   if( ! obj->get_function() ) // the FRealObjective has no Function
    f_lf = nullptr;
   else {
    f_lf = dynamic_cast<LinearFunction *>( obj->get_function() );
    if( ! f_lf )
     throw( std::logic_error( "the objective is not a LinearFunction" ) );

    if( ! f_lf->get_num_active_var() )  // the LinarFunction has no Variable
     f_lf = nullptr;
    }
   }

  v_c05f.resize( sb.size() );

  for( Index i = 0 ; i < sb.size() ; ++i ) {  // for each sub-block

   // the objective function of each sub-block must be a C05Function - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   auto obj = dynamic_cast< FRealObjective * >( sb[ i ]->get_objective() );
   if( ! obj )
    throw( std::logic_error( "the objective is not a real function" ) );

   auto c05f = dynamic_cast<C05Function *>( (obj)->get_function() );
   if( ! c05f )
    throw( std::logic_error( "the objective is not a C05Function" ) );
   v_c05f[ i ] = c05f;

   // nephews are not allowed- - - - - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( sb[ i ]->get_nested_Blocks().size() )
    throw( std::logic_error( "nephew are not allowed" ) );

   // Variable of sub-Block are not expected, neither are Constraint - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( sb[ i ]->get_static_variables().size() )
    throw( std::logic_error( "static Variable are not allowed" ) );

   if( sb[ i ]->get_dynamic_variables().size() )
    throw( std::logic_error( "dynamic Variable are not allowed" ) );

   if( sb[ i ]->get_static_constraints().size() )
    throw( std::logic_error( "static Constraint are not allowed" ) );

   if( sb[ i ]->get_dynamic_constraints().size() )
    throw( std::logic_error( "dynamic Constraint are not allowed" ) );
   }
  }

 // the set of "active" Variable in all Function must be the same- - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NumVar = v_c05f[ 0 ]->get_num_active_var();
 LamVcblr.resize( NumVar );
 auto vi = std::as_const( v_c05f[ 0 ] )->begin();
 for( Index i = 0 ; i < NumVar ; ++vi )
  LamVcblr[ i++ ] = static_cast< ColVariable * >( & (*vi) );

 if( f_lf ) {
  if( f_lf->get_num_active_var() != NumVar )
   throw( std::logic_error( "the list of active Variable do not match" ) );

  auto v = f_lf->begin();
  for( auto vi = LamVcblr.begin() ; vi != LamVcblr.end() ; ++v , ++vi )
   if( static_cast< ColVariable * >( & (*v) ) != *vi )
    throw( std::logic_error( "the list of active Variable do not match" ) );
  }

 for( Index i = 1 ; i < sb.size() ; ++i ) {
  if( v_c05f[ i ]->get_num_active_var() != NumVar )
   throw( std::logic_error( "the list of active Variable do not match" ) );

  auto v = v_c05f[ i ]->begin();
  for( auto vi = LamVcblr.begin() ; vi != LamVcblr.end() ; ++v , ++vi )
   if( static_cast< ColVariable * >( & (*v) ) != *vi ) 
    throw( std::logic_error( "the list of active Variable do not match" ) );
  }
  
  
  Index NumBVar = 0;  // count the number of Variable in the Block
 auto v_s_Variable = f_Block->get_static_variables();
 for( auto & el : v_s_Variable ) {
  if( un_any_thing_0( ColVariable , el , ++NumBVar ) )
   continue;
  if( un_any_thing_1( ColVariable , el , NumBVar += var.size() ) )
   continue;
  if( un_any_thing_K( ColVariable , el , NumBVar += var.num_elements() ) )
   continue;
  throw( std::logic_error( "some static Variable is not a ColVariable" ) );
  }

 if( NumBVar < NumVar )
  throw( std::logic_error( "too few ColVariable in the Block" ) );

 // check that the Variable in the Block agree with that in the C05Function- -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::vector<ColVariable *> LamBVcblr( NumBVar );

 Index cnt = 0;
 for( auto & el : v_s_Variable )
  un_any_static( el , [ & ]( ColVariable & sv ) { LamBVcblr[ cnt++ ] = & sv;
                  } , un_any_type<ColVariable>() );

 std::sort( LamBVcblr.begin() , LamBVcblr.end() );

 std::vector<ColVariable *> LamVcblrO( LamVcblr );
 std::sort( LamVcblrO.begin() , LamVcblrO.end() );

 if( ! std::includes( LamBVcblr.begin() , LamBVcblr.end() ,
		      LamVcblrO.begin() , LamVcblrO.end() ) )
 throw( std::logic_error(
		   "some ColVariable in C05Function are not in the Block" ) );

 LamVcblrO.clear();
 LamBVcblr.clear();

 // no dynaimic variables are allowed  - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_Block->get_dynamic_variables().size() )
  throw( std::logic_error( "dynamic Variable are not allowed" ) );

 // if some Constraint are present, their functions must be linear
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_Block->get_static_constraints().size() )
  throw( std::logic_error( "static Constraint are not allowed" ) );

 if( f_Block->get_dynamic_constraints().size() )
  throw( std::logic_error( "dynamic Constraint are not allowed" ) );
  
  // read information about the function  - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NrFi = v_c05f.size();
 
  // allocate memory- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 Lambda.resize( NumVar );    // the default starting point=0
 LambdaBar.resize( NumVar);// stability center default=0
 
 for(Index i=0; i<NumVar; i++)
	Lambda[i]=LambdaBar[i]=0;
 
 LmbdBst.resize( NumVar );

 Result = ThinComputeInterface::kError;
 DoSS = false;
 
 Lowerbound=-Inf<VarValue>();// lower bound on the function
 TrueLB=false;
 
 FiBest = FiLambda = Inf<VarValue>(); // Fi( Lambda ) is not known
 
 // allocate memory for the subgrads and the direction - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // The direction will be allocated after either Deflection will have been
  // defined or changed [see InitDeflection()]. Instead, the subgradient is
  // considered separately from the deflection and it will keep on having the
  // same values after any deflection update.
  gi=new VarValue[NumVar];
  Gi.resize(NumVar);
  dir.resize(NumVar);

  NrmGi = 0;
  NrmDir=0;
 
 // asking the projection of both DirM1 and and Gi implies the automatic
  // projection of dir- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( SGPar1 & PrjDIRM1 ) && ( SGPar1 & PrjSG ) )
   SGPar1 &= ~PrjDIR;
   
  // reset algorithm  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
  // checking for 0-th component - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index * SG0Bse;
  Index SG0BDm =GetGi( gi , SG0Bse ,0, NumVar);
    for(Index i=0; i<NumVar; i++)
      Gi[i]=gi[i];
  ZeroComp = SG0BDm > 0;
  // initialization of data structures for incremental done there inside

  // initialize HatLambda  - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( SGPar4 ) {
   LambdaHat.resize(NumVar);
   for(Index i=0; i<NumVar; i++)
	LambdaHat[i]=VarValue( 0 );
   }

  FiHat = 0;     // Fi( HatLambda ) is not known

  // initialize both stepsize and deflection solvers - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if(stepsize)
		InitStepsize();
	
	InitDeflection();  

  // other initializations- - - - - - - - - - - - - - - - - - - - - - - - - -
  fs = ThinComputeInterface::kUnEval;
 
  EmptySet = false;       // the feasible set is not empty

/*-----------------------------------------------------------------------------*/
	for(Index i=0; i<NrFi; i++){
		IsConvex = IsConvex && (v_c05f[i]->is_convex()); //is true iff every component of the objective function is convex
		IsConcave = IsConcave && (v_c05f[i]->is_concave()); //is true iff every component of the objective function is concave
	}
	
	UpdtLowerBound();
	
	AllocateUBLB(0, NumVar);

 ReSetAlg( RstAlgPar );  // Lambda is reset inside
 
}//end SubGradSolver::Set_Block

 /*--------------------------------------------------------------------------*/

int SubGradSolver::compute(bool changedvar)
{
 if( MaxIter == 0 )  // no iteration must be performed
  return( kStopIter );
  

 // basic sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! f_Block )
  return( kBlockLocked );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - 
while( ! v_mod.empty() ) {
  bool owned = f_Block->is_owned_by( f_id );       // check if already locked  
  if( ( ! owned ) && ( ! f_Block->read_lock() ) )  // if not try to read_lock
   return( kBlockLocked );                         // return error on failure
   
if( ! owned )             // if the Block was actually read_locked
   f_Block->read_unlock();  // read_unlock it
  }
  
  //Check if the block can be used
  
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// start timer - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 c_start = std::clock();
  
	// initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 double lastETT = 0;  // last "time" eEveryTTime events have been called
 ParIter=0; //number of iterations
 Result = Solver::kOK; //risultato inizializzato a default
 SCalls++; //number of calls to compute() included this one
 Index comp;	//a component

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 do {
  // function log: print the information relative to the stability center
  // and the best point - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log1();

  // check for running time - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( elapsed() > MaxTime ) 
  {
	VLOG(1," ~ stop due to max time" << std::endl );
    Result = kStopTime;
   break;
   }
   
   // run time-periodic events - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( EveryTTm && ( elapsed() >= lastETT + EveryTTm ) ) {
   for( auto & ev : v_events[ eEveryTTime ] ) {
    auto res = ev();
    if( res == eStopOK ) { Result = kOK; break; }
    if( res == eStopError ) { Result = kError; break; }
    }

   lastETT = elapsed();  // reset counter
   }
  
  if( Result != kUnEval ) {
   VLOG( 1 , " ~ stop due to time-periodic event" << std::endl );
   break;
   }
  // To compute the next point Lambda1, select the direction Dir and the
  // stepsize Nu using, respectively, Deflection and Stepsize objects
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  if(FiBar <Inf<VarValue>()){
	 
	//if the function value (at Lambdabar) is unknown  skip the direction and stepsize
   // computation and go straight to FiAndGi()

   // direction log: print the information relative to DIRM1  - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Log2();

   // FormD call: a new direction is provided, but Dir is not updated
   // yet because AddVariables() and RemoveVariables() need to work
   // with the prior vector Dir - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   FormD( );

	// check the status and take the action, only during
   // outer iterations  - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( !InnIter )
    fs = f_lf->compute(false);

   if( fs == ThinComputeInterface::kError ) {
    VLOG( 1 , " ~ STOP" << std::endl );
    Result = Solver::kError;
    break;
    }
	
   // if something is changed go back to FormD(): the deflection coefficient
   // Alfa, and in turn Dir, must be again computed, and, in the
   // deflection-restricted approach, even the stepsize Nu- - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	
   // save the search direction Dir (and calulate the stepsize in the stepsize-restriction scheme)
   
   //questa parte di sotto riguarda qualcosa sul slavataggio in memoria che non ho capito
   
   //: note that the previous DIR is not saved
   // as DirM1, but it is merely discharged; only and one vector is kept
   // for storing the direction- - - - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
   SaveDir();
   
   // if the Stepsize object needs the scalar product between DirM1 and Gi (<d_{i-1}, g_i>),
   // then the direction vector is kept- - - - - - - - - - - - - - - - - - - -

   if( stepsize->NeedsdkM1Gk() && ( ! deflection ) ){
	for(Index i=0; i<NumVar; i++)
		dir[i]=Gi[i];
   }
   // direction log: print the information relative to DIR   - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Log2();
	
	if( step <= 1e-8 * tStar ) { // it'a a shortstep
    CSmallStep++; // another shortstep 
    if( CSmallStep > ( ( NItIncr > 0 )? 100 * NrFi : 100 ) ) { // If the number of incremental iteration after an outer iteration (NItIncr) is >0 then I compare CSmallStep with 100*NrFi (numbers of component of the obiective sum function), else with 100
	// If the number of consecutive shortstep is higher than the limit value, then I stop the iteration and return the status
     Result = kStopIter;
     break;
     }
    }
   else// it's a normal step 
    CSmallStep = 0;//The counter is reinizialized because we only check about consecutive short steps
   }  // end( search direction computation, when there are no info about the function valuated in LambdaBar) - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
  // check for optimality  - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( ! InnIter ) && IsOptimal() ){
	   // run optimality events
   int res = eContinue;
   for( auto & ev : v_events[ eBeforeTermination ] )
    if( ( res = ev() ) != eContinue )
     break;

   if( res == eStopError ) {
    //VLOG( 1 , " ~ stop (error) by optimality event" << std::endl );
    Result = kError;
    break;
    }
    
	if(res==eForceContinue){
		throw("This is not allowed");
		break;
	}
   //VLOG( 1 , " ~ stop (optimal)" << std::endl );
   Result = kOK;
   break;
   }

  // run iteration-periodic events- - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( EverykIt && ( ! ( ParIter % EverykIt ) ) )
   for( auto & ev : v_events[ eEverykIteration ] ) {
    auto res = ev();
    if( res == eStopOK ) { Result = kOK; break; }
    if( res == eStopError ) { Result = kError; break; }
    }

  if( Result != kUnEval ) {
   //VLOG( 1 , " ~ stop due to time-periodic event" << std::endl );
   break;
   }

  // calculate the tentative point Lambda1 (not projected yet)   - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  FormLambda1();

  //VLOG( 1 , std::endl << "           " << " |Lambda1| = " << Norm( Lambda , NumVar )  );

  // a real iteration (iterations where Fi() is not evaluated count as well) -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ParIter++; //new iteration

  // is an inner iteration?- - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  InnIter = NItIncr? GiEvaltns % ( NItIncr + 1 ): 0; 
  //InnIter takes GiEvaltns%(NItIncr+1) if NItIncr!=0, else 0
  //GiEvaltns is the number of calls to Gi()
  
  LHasProj = false; //lambda has not been projected yet in this iteration
  
  // pick up component comp (defined at the beginning of the function) - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if (InnIter){
		while (!comp)// I don't want the first component
			comp=myrandom(NrFi);//I take one of the component
	}
	else
		comp=Inf<Index>();

	// project Lambda1 onto the feasible set - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	while (1){
		if(ProjFsb(Lambda, EmptySet)){//check if lambda components were changed in the projection
			LHasProj=LHasChgd=true;
			//VLOG( 1 , " ~ |Lambda1| = " << Norm( Lambda , NumVar )  );
			if(EmptySet)
				break;
			//continue while the set is feasible: so it will be possible to take a point	
			}
	
		// never evaluate the objective function at the same point - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		if( ! LHasChgd )//if the point lambda is the same as the previous iteration
			throw ("SubGradSolver::compute() this should not happen");
	
		// compute Fi( Lambda1 ) and a subgradient / constraint Gi[wFi]- - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		FiAndGi( comp );//update LmbdBst
   
		// if the Lambda1 is infeasible, go back to the projection - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		if( FiLambda < Inf<VarValue>() )
			break;

   			} // end function evaluation - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	
	// project the subgradient Gi[wFi] over the tangent cone at tentative point
   // Lambda1- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( SGPar1 & PrjSG )// check if Gi[wFi] is not projected yet
   ProjTng( Gi , Lambda );

  NrmGi = Norm( Gi , NumVar );//compute the sqared norm of Gi[wFi]

  // compute the scalar product between Dir and Gi[wFi] at Lambda1- - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
  if( stepsize->NeedsdkM1Gk() ) {
   dM1Gk = ScalarProduct( dir , Gi , NumVar );
   dM1GkDone = true;
   }
  else
   dM1GkDone = false;
   
    // print the logs regarding Fi and Gi[wFi] at Lambda1 - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   //log su fi(lambda1) e gi(wfi) 
   
   // check for feasibility- - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
   if( EmptySet ) {
		Result = kInfeasible;
   break;
   }
   
   if (FiLambda==-Inf<VarValue>()){
	   Result=kUnbounded;
	   break;
   }
  // check whether or not the Lower Bound has changed - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
	UpdtLowerBound();
  
  //VLOG2( 1 , LowerBound > - Inf<VarValue>() , " ~ UB = " << -LowerBound );
   
   // if FiBar (the function at LambdaBar) is unknown, or no deflection is
  // used, then move to Lambda1 - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( ( FiBar == Inf<VarValue>() ) || ( ! deflection ) ) {
		
		GotoLambda1();
		DoSS = true;//A serious step is done
		//VLOG( 1 , std::endl );
		continue;
   }
   
   // Null Step (NS) vs Serious Step (SS)- - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  DoSS = deflection->DoSS();
  if( DoSS ) {  // SS- - - - - - - - - - - - - - - -

   VLOG2( 1 , deflection->Delta() < Inf<VarValue>() , std::endl << " SS[" << CSSCntr );

   GotoLambda1();

   CSSCntr++;//new serious step is done
  CNSCntr = 0;// turn the number of consecutive null steps to 0
  }else {  // NS- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   VLOG2( 1 , deflection->Delta() > - Inf<VarValue>() , std::endl << " NS[" << CNSCntr );

   UpdateSigma();//need to modify the error 

   CNSCntr++;
   CSSCntr = 0;

   }  // end SS/NS - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VLOG( 1 , std::endl );

  // check for Lower Bound stopping condition - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( TrueLB )
   if( FiBest - EpsLin * std::abs( FiBest ) <= Lowerbound )
    break;

  if( FiBest <= Lowerbound  * ( 1.0 - ( Lowerbound > 0.0 ? EpsLin : - EpsLin )) ) {
   Result = kUnbounded;
   break;
   }
   
  } while( ( ! MaxIter ) || ( ParIter < MaxIter ) );
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  return( Result );

}//end ( SubGradSolver::compute() )

/*--------------------------------------------------------------------------*/ 
 
 /*--------------------------------------------------------------------------*/
 
void SubGradSolver::set_par( const idx_type par , const int value ){
	if(par<FirstintDeflectionPar){
 switch( par ) {
  case( CDASolver::intMaxIter ):
   if( value < 0 )
    throw( std::invalid_argument( "MaxIter must be >= 0" ) );
   MaxIter = value;
   break;
  case( CDASolver::intLogVerb ):
   LogVerb = value;
   break;
   case( Solver::intEverykIt ):
	EverykIt=value;
  case(intSGPar1):
		if(value<1)
		throw( std::invalid_argument( "SGPar1 must be >= 1" ) );
		SGPar1 = value;
	break;
   case(intSGPar3):
		if(value<1)
		throw( std::invalid_argument( "SGPar3 must be >= 1" ) );
		SGPar3 = value;
	break;
	case(intSGPar4):
		SGPar4=value;
	break;
	case (intSGPar5):
		SGPar5=value;
	break;
	case( intRstAlgPar ):
		RstAlgPar = value;
  break;
	default:
		CDASolver::set_par( par , value );
	}
}
else{
	 if(par<FirstintStepsizePar)
		deflection->set_par_d(par-FirstintDeflectionPar, value);
	 else
		stepsize->set_par(par-FirstintStepsizePar,value);
	}
 }  // end( SubGradSolver::set_par( int ) )
 
 /*--------------------------------------------------------------------------*/
  
void SubGradSolver::set_par( const idx_type par , const double value ){
	if(par<FirstdblDeflectionPar){
		switch( par ) {
			case( dblSGPar2 ):
				if( value < 1 )
					throw( std::invalid_argument( "SGPar2 must be >= 1" ) );
				SGPar2 = value; 
			break;
			case (dbltStar):
				tStar=value;
			break;
			case(dblEDcrs):
				EDcrs=value;
			break;
			case(dblEFnalt):
				EFnal=value;
			break;
			case(dblEInit):
				EInit=value;
			break;
			case (dblEStps):
				EStps=value;
			break;
			case(dblEpsLin):
				EpsLin=value;
			break;
			case( dblMaxTime ):
				if( value <= 0 )
				throw( std::invalid_argument( "dblMaxTime must be > 0" ) );
				MaxTime = value;
			break;
			default:
				CDASolver::set_par( par , value );
		}
	}
	else{
		if(par<FirstdblStepsizePar)
			deflection->set_par_d(par-FirstdblDeflectionPar,value);
		else
			stepsize->set_par(par-FirstdblStepsizePar, value);
	}
} // end( SubGradSolver::set_par( double ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void SubGradSolver::set_log( std::ostream *log_stream){
 f_log=log_stream;
  if( LogVerb > 1 ) {
   *f_log << std::endl << "SubGradSolver: Vars = " << NumVar << " ("
	   << ") ~  NItIncr = " << NItIncr <<  " ~ EpsLin = "
	   << EpsLin << " ~ FiEps = " << EInit << " ~ t* = " << tStar;
   *f_log << std::endl << " SGParameters = ( " << SGPar1 << " ~ " << SGPar2
	   << " ~ " << SGPar3 << " ~ " << ( SGPar4? "true" : "false")
	   <<  " ~ "  << SGPar5 << " )" << std::endl;
   }
  }  // end( SubGradSolver::set_log() )

Index SubGradSolver::GetNumVar( void ) {
  	return(NumVar);
  }

/*--------------------------------------------------------------------------*/

Index SubGradSolver::GetNrFi( void ){
   	return(NrFi);
   	}

/*--------------------------------------------------------------------------*/

int SubGradSolver::get_int_par( const idx_type par ) {
	if(par<intLastParSG){
 switch( par ) {
	case( intMaxIter ):
		return( MaxIter );
	break;
	case( intLogVerb ):
		return( LogVerb );
	break;
	case( intEverykIt ):
		return( EverykIt );
	break;
	case(intSGPar1):
		return(SGPar1);
	break;
	case(intSGPar3):
		return(SGPar3);
	break;
	case(intSGPar4):
		return(SGPar4);
	break;
	case(intSGPar5):
		return(SGPar5);
	break;
	case( intRstAlgPar ):
		return(RstAlgPar);
    break;
	default:
		return(CDASolver::get_dflt_int_par(par));
	}}
	else{
		if(par<FirstintDeflectionPar)
		return(deflection->get_int_par_d(par-FirstintDeflectionPar));
		else
		return(stepsize->get_int_par(par-FirstintStepsizePar));
	}
 }//end (SubGradSolver::get_int_par)

/*--------------------------------------------------------------------------*/

double SubGradSolver::get_dbl_par( const idx_type par ) {
if(par<dblLastParSG){
	switch( par ) {
		case( dblMaxTime ):
			return( MaxTime );
		break;
		case(dblSGPar2):
			return(SGPar2);
		break;
		case(dbltStar):
			return(tStar);
		break;
		case(dblEDcrs):
			return(EDcrs);
		break;
		case(dblEFnalt):
			return(EFnal);
		break;
		case(dblEInit):
			return(EInit);
		break;
		case(dblEStps):
			return(EStps);
		break;
		case(dblEpsLin):
			return(EpsLin);
		break;
		default:
			return( CDASolver::get_dflt_dbl_par( par ) );
		}	
	}else{
		if(par<FirstdblStepsizePar)
			return(deflection->get_dbl_par_d(par-FirstdblDeflectionPar));
		else
			return(stepsize->get_dbl_par(par-FirstdblStepsizePar));
	}
}//end (SubGradSolver::get_dbl_par)

void SubGradSolver::ReSetAlg( unsigned char RstLvl){
 if( ! ( RstLvl & RstAlg ) ) {  // reset algorithmic parameters - - - - - - -
  ParIter = 0;             // reset iterations counerst, comprised
  CSmallStep = 0;          // consecutive NS/SS count
  CSSCntr = CNSCntr = 0;
  EInit=dflt_dbl_par[dblEInit-dblLastParCDAS];  // reset Fi precision
  }

 if( ! ( RstLvl & RstCrr ) ) {
	// Vec_VarValue z[NumVar];
	  SetLambda(Vec_VarValue()); 
 } // reset the current point to all-0- - - - - - -


 if( ! ( RstLvl & ( RstSbg ) ) ) {  // reset everything - - - - - - -
  guts_of_destructor();
  int x=0;
  }
 else
  if( ! ( RstLvl & RstSbg ) ) {  // reset the subgradient and the directions
	                         // (but not the knapsack constraints)- - - -
   //for( Index i = MaxName ; i-- > MaxNCnst ; )
    //Oracle->Deleted( i );
    int j=0;
   }

 if( ! ( RstLvl & RstFiV ) )  // reset Fi( LambdaBar )- - - - - - - - - - - -
  FiBar = Inf<VarValue>();

 }  // end( SubGradSolver::ReSetAlg() )
 
 
/*--------------------------------------------------------------------------*/
/*-----------------------  METHODS FOR READING THE RESULTS --------------------------*/
/*--------------------------------------------------------------------------*/
Vec_VarValue SubGradSolver::ReadBestSol( void ){
 //I = nullptr;
 //D = NumVar;
 return( KpBstL ? LmbdBst : LambdaBar );

 }  // end( SubGradSolver::ReadBestSol() )

/*--------------------------------------------------------------------------*/

VarValue SubGradSolver::ReadFiVal( void ){
 //if( wFi!= Inf<Index>() )
  //throw( std::invalid_argument( "ReadFiVal: wFi != INF" ) );
    return(FiBar);
 }  // end( SubGradSolver::ReadFiVal() )

/*--------------------------------------------------------------------------*/

Vec_VarValue SubGradSolver::ReadSol( void ){
 return( LambdaBar );

 }  // end( SubGradSolver::ReadSol() )

/*--------------------------------------------------------------------------*/

VarValue SubGradSolver::ReadBestFiVal( void ){
 //if( wFi!= Inf<Index>() )
  //throw( std::invalid_argument( "ReadBestFiVal: wFi != INF" ) );
	if(IsConvex)
 		return( FiBest );
 	else {
 		if(IsConcave)
 			return (-FiBest);
 		else
 			throw (std::invalid_argument("ReadBestFiVal: Fi() is neigther convave or convex"));
 		}

 }  // end( SubGradSolver::ReadBestFiVal() )

/*--------------------------------------------------------------------------*/

VarValue SubGradSolver::ReadHatFiVal( void ){
 if( SGPar4 ){
 	if(IsConvex)
 		return( FiHat );
 	else {
 		if(IsConcave)
 			return (-FiHat);
 		else
 			throw (std::invalid_argument("ReadBestFiVal: Fi() is neigther convave or convex"));
 		}
 }
 else{
 	if(IsConvex)
 		return( Inf<VarValue>() );
 	else {
 		if(IsConcave)
 			return (-Inf<VarValue>());
 		else
 			throw (std::invalid_argument("ReadBestFiVal: Fi() is neigther convave or convex"));
 		}
 
 }
 }  // end( SubGradSolver::ReadHatFiVal() )


bool SubGradSolver::IsOptimal( VarValue eps ){
 VarValue FiL = FiBest;

 if( FiL == Inf<VarValue>() )
  return( false );
 else {
  if( FiL < 0 ) FiL = - FiL;
  if( FiL < 1 ) FiL = 1;

  if( eps <= 0 )
   eps = EpsLin;

  if( SGPar4 )
   return( tStar * std::sqrt( NrmDir ) + std::min( Epsilon , HatEpsilon )
	   <= eps * FiL );
  else
   return( tStar * std::sqrt( NrmDir ) + Epsilon <= eps * FiL );
  }
 }  // end( SubGradSolver::IsOptimal() )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

 void SubGradSolver::SetLambda(Vec_VarValue tLambda){
	//assign tLambda to Lambdabar - - - - - - - - - - - - - - - - - - - - -
	Index i=0;
	if(tLambda.size()>0){
		for(i=0; i<NumVar; i++)
			LambdaBar[i]=tLambda[i];
	}
	else{
		for(i=0;i<NumVar; i++)
			LambdaBar[i]=VarValue(0);
	}
	FiBest=FiBar=Inf<VarValue>();
	
}// end (SubGradSolver::SetLambda())

/*--------------------------------------------------------------------------*/

void SubGradSolver::FormD(void){
	//d_i=a_i*G_i+(1-a_i)D_{i-1}
	//a_i stepsize
	//G_i=P_{T_I}(g_i) the projection of the subgradient over the tangent cone 
	//D_{i-1}=P_{T_i}(d_{i-1}) the projection over the tangent cone of the previous direction
	//d_{-1}=0
/*--------------------------------------------------------------------------*/
/*Check if the component are projected over the tangent cone according to the SGPar1 parameter:
it's coded bitwise and every bit is equal to 1 if the field is projected, 0 otherwise
-bit 1 -> g_i
-bit2 ->d_{i-1}
-bit 3 ->d_i
*/
if ((SGPar1 & PrjSG) && (! DoSS) ) //status check(fs==//status check)) //g_i is already projected but no serious step has been done
//The projection has to be recomputed
{
	ProjTng(Gi, LambdaBar);
NrmGi=Norm( Gi, NumVar);// the squared norm
dM1GkDone = false;//the scalar product \f$d_{i-1}^{\top} g_i\f$ has to be computed
  }
// subgradient log: print the information relative to Gi[wFi] - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Log4();

 // project DirM1 and compute its norm - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ((SGPar1 & PrjDIRM1)) //status check(fs==//status check))
	//if d_{i-1} has been projected
	{if(SGPar1 & PrjDIR){
		//if d_i has been projected over the tangent cone of the previous Lambdabar
		if (DoSS){ //if a serious step has been done
			ProjTng(dir, LambdaBar);// mandatory to reproject
			NrmDir=Norm(dir, NumVar);
			dM1GkDone = false; //the scalar product \f$d_{i-1}^{\top} g_i\f$ has to be computed
		}	
	}
	else{//is needed to project d_i
		ProjTng(dir, LambdaBar);
		NrmDir=Norm(dir, NumVar);
		dM1GkDone = false; //the scalar product \f$d_{i-1}^{\top} g_i\f$ has to be computed
	}
}	

// (re)compute the scalar product between DirM1 and Gi[wFi] at Lambda - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (deflection){
	if(! dM1GkDone)//if the scalar product has not been already computed
		dM1Gk=ScalarProduct( Gi, dir, NumVar);
}
else{
	//no deflecton is done -> d_i=g_i
	dGk=NrmDir=NrmGi;// the scalar product <d_i, g_i>
	if(! dM1GkDone)
dM1Gk=dGk;// <d_{i-1}, g_i>=<d_i, g_i>
}

dM1GkDone=true; //the scalar product <d_{i-1}, g_i> has been done 

/*-------------------------------------------------------------------------------*/
/* Check on the approach: stepsize-restricted or deflection-restricted.
This information is coded bit wise in SGPar3: second bit is equal to 1
in the first approach, equal to 0 in the secon one*/
if (SGPar3 & StepRST){
	//stepsize-restricted approach: deflection, compute, stepsize
	if (deflection){
		//deflection computation
		deflection->NewDEF();// inizialization
		alpha=deflection->GetDFLCoeff();// get alpha 
	}
	//vlog su chi è alpha
}
else{
	//deflection-restricted approach: stepsize, deflection, compute
	dGk=dM1Gk;// is not important the distinction between them
	
	if( ! InnIter)// if is an outer iteration is needed to recompute stepsize
	stepsize->NewStep();
	
	step=stepsize->GetStepsize(InnIter);// get the stepsize according to the type of iteration
	
	//vlog sul valore di step
	
	deflection->NewDEF();// the case of NO DEFLECTION is not allowed
	alpha=deflection->GetDFLCoeff();// get alpha
	
	//vlog sul valore di alpha
	/*-------------------------------------------------------------------------------------------*/
	
	/* Check if any safe rule is taken according to SGPar3 coded bit-wise: 
	the first bit is equal to 1 if the rule is used, it's equal to 0 otherwise*/
	
	if((SGPar3 & SafeRULE) && !InnIter){
		//if a safe rule is used and there is an outer iteration
		VarValue target= (step * NrmDir)/(FiBar-stepsize->GetLev()+ (step + NrmDir)); //(v_i*||d_{i-1}||^2)/(f_i-f_i^lev + v_i*||d_{i-1}||^2 )
		//f_i^lev=f_i^rec-\delta_i
		//f_i^rec==min{f_k:k=1,..,1}
		//\delta_i>0 depends on the choice of vanish or nonvanishing method
		//perchè diviso e non *?????
		alpha=std::min(1.0, std::max(alpha, target));// to respect the right inequality for the safe rule
		
		//vlog sul valore dell'alpha appena calcolato
	}
}

// print Dir norm - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 //Capire cosa ci va eventualmente

}// end SubGradSolver::FormD()

/*--------------------------------------------------------------------------*/
void SubGradSolver::SaveDir(void){
// compute the linearization error of Dir at LambdaBar- - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if(alpha==1)
	Epsilon=Sigma;
 else
	Epsilon= alpha* Sigma + (VarValue(1)-alpha)*Epsilon;
 
  // change LambdaHat - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( SGPar4 )// is a boolean type: it controls if LambdaHat is taken
  ChgLambdaHat();

 // save the search direction Dir  - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if (deflection){
	 //d_i is the convex combination of d_{i-1} and g_i
	 
	 if(alpha==1){
		 for(Index i=0; i<NumVar; i++)
			dir[i]=Gi[i];
	 }
	 else{
		 for(Index i=0; i<NumVar; i++){
			 VarValue g=dir[i];//if we do the k-th iteration this is d_{k-1}[i]
			 dir[i]=alpha*Gi[i]+(VarValue(1)-alpha)*g;
		 }
			
	 }
	 
	 //Ci sarebbe lo swap, delle posizioni, vedi fine FormD()
	 
	 // project Dir on tangent cone at LambdaBar  - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if(SGPar1 & PrjDIR)// capire perchè fa questo controllo
		ProjTng(dir, LambdaBar);
	
	// compute the (squared) norm of the direction and its scalar product with
  // the subgradient   - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	NrmDir = Norm( dir , NumVar );
	dGk = ScalarProduct( Gi , dir , NumVar );

  }  // end( deflection )  - - - - - - - - - - - - - - - - - - - - - - - - - -
	
	// stepsize computation - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( SGPar3 & StepRST ) {//check on the approach: stepsize restricted
  if( ! InnIter ) // newstep only in outer iterations
   stepsize->NewStep( );

  step = stepsize->GetStepsize( InnIter );
  
  // VLOG( 1 , " ~ step = " << step );
  // safe rule - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( SGPar3 & SafeRULE ) {//check on the safe rule
   stepsize->SetMaxBeta( alpha );//capire
   step = stepsize->GetStepsize( InnIter );
   // VLOG( 1 , " ~ (fixed) step = " << step );
   
   }

  }  // end stepsize-restricted scheme - - - - - - - - - - - - - - - - - - - -
  
 }  // end( SubGradSolver ::SaveDir() )

/*--------------------------------------------------------------------------*/
void SubGradSolver::FormLambda1(void){
	if(deflection){
		for( Index i = 0 ; i < NumVar ; i++ )
			Lambda[i]=LambdaBar[i] - step * dir[i];
	}
	else{
		for( Index i = 0 ; i < NumVar ; i++ )
			Lambda[ i ] = LambdaBar[ i ] - step * Gi[ i ];
	}
 	LHasChgd = true;// lambda has changed since the last call to FiandGi()
}	//end(SubGradSolver::FormLambda1())

/*--------------------------------------------------------------------------*/
void SubGradSolver::FiAndGi(  Index wFi ){
 // pass Lambda1 - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( LHasChgd ) {  // if Lambda1 has changed, set the new one 
  SetLambda( Lambda );
  LHasChgd = false;
  }

 // if outer iteration, evaluate Fi in Lambda1  - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( wFi == Inf<Index>() ) {  // call Fi()- - - - - - - - - - - - - - - - - -
  FiLambda = f_lf->get_value();//f(lambda)
  GiEvaltns++;

  if( FiLambda == -Inf<VarValue>() )    // Fi() unbounded (-inf)
   return;

  // update FiBest, if necessary - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( FiLambda < FiBest ) {
   FiBest = FiLambda;
   if( KpBstL )
    for(Index i=0; i<NumVar; i++)
		LmbdBst[i]=Lambda[i];
   }
  }  // end( function evaluation) - - - - - - - - - - - - - - - - - - - - - - -

 // if Lambda is unfeasible, a new subgradient is computed; otherwise a new
 // constraint is added - - - - - - - - - - - - - - - - - - - - - - - - - - -

if( FiLambda == Inf<VarValue>() )// it is a constraint
	throw("It's out of the feasible region");
 else
  EvalGi( wFi );  // the item is a subgradient 

}//end(SubGradSolver::FiandGi())

/*--------------------------------------------------------------------------*/
void SubGradSolver::GotoLambda1( void ){
 // compute the linearization error of Dir at the new LambdaBar  - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiBar < Inf<VarValue>() && deflection ) {
  VarValue SclPrd;
  if( ! LHasProj )
   SclPrd = step * NrmDir;
  else {
   SclPrd = 0;
   for( Index i = 0 ; i < NumVar ; i++ )
    SclPrd += dir[ i ] * ( LambdaBar[ i ] - Lambda[ i ] );
   }

  Epsilon += FiLambda - FiBar + SclPrd;
  }

 // move Lambda to LambdaBar= - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	for(Index i=0; i <NumVar; i++)
		LambdaBar[i]=Lambda[i];
	if( ! InnIter )
		FiBar = FiLambda;

 }  // end( SubGradSolver::GotoLambda1() )

/*--------------------------------------------------------------------------*/
Index SubGradSolver::GetGi(VarValue *sub, Index * sgb, Index strt, Index stp){
	Index NumVar = v_c05f[ 0 ]->get_num_active_var();
	if( stp > NumVar )
		stp = NumVar;
		auto range=std::make_pair(strt, stp);
	HasLin=f_lf->has_linearization(true);
	if(HasLin){
		f_lf->get_linearization_coefficients(sub, range);
		sgb=nullptr;
	}
	else
		throw("It's impossible to compute the subgradient because only diagonal linearization is allowed");

	return(stp-strt);
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void SubGradSolver::guts_of_destructor( void )
{
  LmbdBst.clear();
  LambdaBar.clear();
  Lambda.clear();
  if(SGPar4)
	  LambdaHat.clear();
  
  LamVcblr.clear();
  Lbv.clear();
  Ubv.clear();
  v_c05f.clear();
 }  // end( SubGradSolver:guts_of_destructor )
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
void SubGradSolver::Log1(void){
	if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;
  
   *f_log << std::endl << "{" << SCalls << "-" << ParIter << "} ";

   *f_log << " FiC = ";
    if( FiBar == Inf<VarValue>() )
     *f_log << " -INF";
    else
      *f_log << -FiBar;

    *f_log << " ~ FiBest = ";
    if( FiBest == Inf<VarValue>() )
     *f_log << " -INF";
    else {
     *f_log << -FiBest;

     if( Lowerbound > - Inf<VarValue>() )
      *f_log << " ~ Gap = " << ( FiBest - Lowerbound ) /
                                std::max( std::abs(Lowerbound ) , 1.0 ) << " %";
     }
}// (end SubGradSolver::Log1())

/*--------------------------------------------------------------------------*/

void SubGradSolver::Log2(void){
	// information relative to the direction
	if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;

 *f_log << std::endl << "            ";
 
 if( SGPar4 )
    *f_log << " ( |D| , Epsilon , HatEpsilon ) = (" <<sqrt(NrmDir)
	    << " , " << Epsilon << " , " << HatEpsilon <<  ") ~ GixD = "
	    << dGk;
   else
    *f_log << " ( |D| , Epsilon ) = (" << sqrt(NrmDir) << " , "
	    << Epsilon <<  ") ~ GixD = " << dGk;

 }  // end( SubGradSolver::Log2() )

/*--------------------------------------------------------------------------*/

void SubGradSolver::Log3( void ){
	if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;

 *f_log << std::endl << "            ";
 
 *f_log << " Fi1 = ";  // print function value
   if( FiLambda == Inf<VarValue>() )
     *f_log << " -INF => STOP.";
   else
    *f_log << -FiLambda;

 }  // end(SubGradSolver::Log3() )

/*--------------------------------------------------------------------------*/

void SubGradSolver::Log4( void )
{
  if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;
   // if the deflection is used the information relative to the subgradient
   // have already been displayed [see Log2() ]

   *f_log << std::endl << "           ";
   if( deflection )
    *f_log << " ( |Gi| , Sigma ) = (" << sqrt(NrmGi) << " , " << Sigma
	    << ")";
 }  // end( SubGradSolver::Log4() )


/*--------------------------------------------------------------------------*/


void SubGradSolver::UpdateSigma( void ){
 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! deflection )
  throw( "UpdateSigma: this is not allowed");

 // compute the linearization error of Gi[wFi] at Lambdabar - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Sigma += FiBar - FiLambda;
 if( ! LHasProj )
  Sigma -= step * ScalarProduct( Gi , dir , NumVar );
 else {
  VarValue SclPrd = 0;
  for( Index i = 0 ; i < NumVar ; i++ )
   SclPrd +=  Gi[ i ] * ( Lambda[ i ] - LambdaBar[ i ] );
  Sigma += SclPrd;
  }

 //VLOG( 1 , " ~ Sigma1 = " << Sigma );

 }  // end( SubGradSolver::UpdateSigma() )

/*--------------------------------------------------------------------------*/

void SubGradSolver::UpdtLowerBound(void){
	//if the objective function is concave i take as lower bound its -upper bound 
	VarValue Lb= IsConvex ? f_Block->get_valid_lower_bound(false) : -(f_Block->get_valid_upper_bound(false));
		if(Lb> -Inf<VarValue>()){
			Lowerbound=Lb;
			TrueLB=true;
		}
	else {
		Lb= IsConvex ? f_Block->get_valid_lower_bound(true) : -(f_Block->get_valid_upper_bound(true));
		Lowerbound=Lb;
		TrueLB=false;
		}		
 }  // end( SubGradSolver::UpdtLowerBound() )
 
/*--------------------------------------------------------------------------*/

 void SubGradSolver::AllocateUBLB( c_Index strt , c_Index num ){
 bool LBc = Lbv.empty();
 bool UBc = Ubv.empty();
 std::vector<VarValue> tUbv;
 std::vector<VarValue> tLbv;
 if(UBc)
	tUbv.resize(num);
 else
	tUbv.resize(Ubv.size()+ strt);
 if(LBc)
	tLbv.resize(num);
 else
	tLbv.resize(Lbv.size()+strt);

 tLbv.assign(num,- Inf<double>());

 Index i = 0;
 if( ( ! LBc ) && ( ! UBc ) )
	  for( ; i < num ; ++i ) {
	   if( static_cast<ColVariable*>(LamVcblr[i])->is_negative()) {
			tLbv[ i ] = 0;
			LBc = true;
	    }

		if( ( tUbv[ i ] = static_cast<ColVariable*>(LamVcblr[i])->get_ub() ) < Inf<double>() )
			UBc = true;

		if( LBc || UBc )
			break;
   }

 if( LBc && ( ! UBc ) )
	  for( ; i < num ; ++i ) {
		   if( static_cast<ColVariable*>(LamVcblr[i])->is_negative())
			  tLbv[ i ] = 0;

			  if( ( tUbv[ i ] = static_cast<ColVariable*>(LamVcblr[i])->get_ub() ) < Inf<double>() )
			UBc = true;
			
			if( UBc )
				break;
   }

 	 if( ( ! LBc ) && UBc )
	  for( ; i < num ; ++i ) {
		  if( static_cast<ColVariable*>(LamVcblr[i])->is_negative()) {
			  tLbv[ i ] = 0;
			  LBc = true;
	    }

	   tUbv[ i ] = static_cast<ColVariable*>(LamVcblr[i])->get_ub();

	   if( LBc )
	    break;
	   }
//
	for( ; i < NumVar ; i++ ) {
		 if( static_cast<ColVariable*>(LamVcblr[i])->is_negative())
			tLbv[ i ] = 0;
		
		tUbv[ i ] = static_cast<ColVariable*>(LamVcblr[i])->get_ub();
	  }

	if( LBc ) {
		if( ( strt == 0 ) && ( num == NumVar ) )
			Lbv = tLbv;
		else {
			Lbv.assign(NumVar+strt,- Inf<double>());
			Lbv.assign(tLbv.begin(),tLbv.begin()+num);

			tLbv.clear();
					}
	}
//
	if( UBc ) {
		if( ( strt == 0 ) && ( num == NumVar ) )
			Ubv = tUbv;
		else {
			Ubv.assign(NumVar+strt,Inf<double>());
			Ubv.assign(NumVar+strt,- Inf<double>());
			Ubv.assign(tUbv.begin(),tUbv.begin()+num);

			tUbv.clear();
		}
	  }
 }  // end( SubGradSolver::AllocateUBLB() )

/*--------------------------------------------------------------------------*/

void SubGradSolver::EvalGi( Index wFi ){
 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( wFi > NrFi )
  wFi = Inf<Index>();//It must be one of the possible index


 Index* SGBse;
 Index SGBDim = GetGi( gi , SGBse , wFi? Inf<Index>() : NumVar, NumVar );
  for(Index i=0; i<NumVar; i++)
    Gi[i]=gi[i];

 // compute the linearization error of Gi[wFi] at Lambda1 - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 VarValue upest=f_lf->get_upper_estimate();
 VarValue lowest=f_lf->get_lower_estimate();
 Sigma = upest-lowest;	
 
 if( SGPar4 )
  SigmaHat = Sigma;


 GiEvaltns++;

 }  // end( SubGradSolver::EvalGi() )
/*--------------------------------------------------------------------------*/
/*------------------------ Projection methods ------------------------------*/
/*--------------------------------------------------------------------------*/

bool SubGradSolver::ProjFsb( Vec_VarValue Lpt , bool & EmptySet ){
 bool LptHasChgd = false;  // it is true, if Lpt changes
 EmptySet = false;         // it is true, if the feasible set is empty
 
  int low=Lbv.size();
  int up=Ubv.size();

 if( ( low ) && ( up ) )
  return( LptHasChgd );

 // project Lpt on the box constraints  - - - - - - - - - - - - - - - - - -

 if( (!low) && (!up) ) {
  // meanwhile, check that they are valid

  Index i = 0;
  for( ; i < NumVar ; i++ ) {
   if( Ubv[ i ] - Lbv[ i ] < EpsLin ) {		
    EmptySet = true;
    break;
    }

   if( Lpt[ i ] > Ubv[ i ] ) {
    Lpt[ i ] = Ubv[ i ];
    LptHasChgd = true;
    break;
    }
   else
    if( Lpt[ i ] < Lbv[ i ] ) {
     Lpt[ i ] = Lbv[ i ];
     LptHasChgd = true;
     break;
     }
   }

  if( ! EmptySet )
   for( ; i < NumVar ; i++ ) {
    if( Ubv[ i ] - Lbv[ i ] < EpsLin ) {	
     EmptySet = true;
     break;
     }

    if( Lpt[ i ] > Ubv[ i ] )
     Lpt[ i ] = Ubv[ i ];
    else
     if( Lpt[ i ] < Lbv[ i ] )
      Lpt[ i ] = Lbv[ i ];
    }
  }
 else
  if( low ) {
   Index i = 0;
   for( ; i < NumVar ; i++ )
    if( Lpt[ i ] < Lbv[ i ] ) {
     Lpt[ i ] = Lbv[ i ];
     LptHasChgd = true;
     }
   for( ; i < NumVar ; i++ )
    if( Lpt[ i ] < Lbv[ i ] )
     Lpt[ i ] = Lbv[ i ];
   }
  else {
   Index i = 0;
   for(  ; i < NumVar ; i++ )
    if( Lpt[ i ] > Ubv[ i ] ) {
     Lpt[ i ] = Ubv[ i ];
     LptHasChgd = true;
     }
   for( ; i < NumVar ; i++ )
    if( Lpt[ i ] > Ubv[ i ] )
     Lpt[ i ] = Ubv[ i ];
   }

 return( LptHasChgd );

 }  // end( SubGradSolver::ProjFsb() )

/*--------------------------------------------------------------------------*/

void SubGradSolver::ProjTng( Vec_VarValue Gpt , Vec_VarValue Lpt , Index strt , Index stp ){
 // project Gpt in the tangent cone at Lpt  - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( stp > NumVar )
  stp = NumVar;
  
  int low=Lbv.size();
  int up=Ubv.size();
	
 if( ( low ) && ( up) )
  return;

 // project Gpt on the box constraints  - - - - - - - - - - - - - - - - - -

 if( (!low) && (!up) ) {
  for( Index j = strt ; j < stp ; j++ )
   if( ( ( Lpt[ j ] >= Ubv[ j ] ) && ( Gpt[ j ] <= 0 ) ) ||
       ( ( Lpt[ j ] <= Lbv[ j ] ) && ( Gpt[ j ] >= 0 ) ) )
    Gpt[ j ] = 0;
  }
 else{
  if( !low ) {
   for( Index j = strt ; j < stp ; j++ )
    if( ( Lpt[ j ] <= Lbv[ j ] ) && ( Gpt[ j ] >= 0 ) )
     Gpt[ j ] = 0;
   }
  else
   for( Index j = strt ; j < stp ; j++ )
    if( ( Lpt[ j ] >= Ubv[ j ] ) && ( Gpt[ j ] <= 0 ) )
     Gpt[ j ] = 0;
 }
 }  // end( SubGradSolver::ProjTng() )

/*--------------------------------------------------------------------------*/
VarValue SubGradSolver::ScalarProduct( Vec_VarValue A, Vec_VarValue B, Index n){
	// returns Sum{ i = 0 .. n - 1 } A[ i ] * B[ i ]
	VarValue s=0; 
	for(Index i=0; i<n; i++)
		s+=(A[i]*B[i]);
	return s;
}

/*--------------------------------------------------------------------------*/

VarValue SubGradSolver::Norm(Vec_VarValue A, Index n){
	VarValue loc=0;
	for(Index i=0; i<n; i++)
		loc=A[i]*A[i];
	return loc;
}

/*--------------------------------------------------------------------------*/

void SubGradSolver::InitStepsize( void ){
 step = 0;
 stepsize=Stepsize::new_Stepsize("FumeroTV", this);
 stepsize->Format();

 }  // end( SubGradSolver::InitStepsize() )

/*--------------------------------------------------------------------------*/

void SubGradSolver::InitDeflection( void ){
 // set to 0 all the information connected with the direction- - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NrmDir = dGk = dM1Gk = 0;
 Epsilon = HatEpsilon = 0;
 alpha = 1;
 dM1GkDone = false;
 Index i=0;

 // initialize the direction dir to zero - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	for(i=0;i<NumVar;i++)
		dir[i]=VarValue(0);

 // if the pointer to Deflection object is not 0, invoke Format() in
 // Deflection, otherwise some conditions have to be taken into account:
 // i) the stepsize-restricted scheme must be used and ii) the projection
 // of the direction is useless
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( deflection ) {
  //DirPos = second;
  deflection->Format();
  }
 else {
  SGPar3 |= StepRST;//or bitwise
  SGPar1 &= ~( PrjDIRM1 | PrjDIR );
  }
 }  // end( SubGradSolver::InitDeflection() )

/*--------------------------------------------------------------------------*/

void SubGradSolver::ChgLambdaHat( void ){
// compute the LambdaHat=lambdahat_{i+1}=alpha_{i+1}*lambda_i+(1-alpha_{i+1})*lambdahat_i- - - - - - - - - - - - - - - -  - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( alpha == 1 ) {
  HatEpsilon = SigmaHat;
  for(Index i=0; i<NumVar; i++)
	LambdaHat[i]=Lambda[i];
  FiHat = FiLambda;//f(lambdahat_i)=f(lambda_i)
  }
 else {
  // compute the linearization error HatEpsilon of DirM1 in LambdaHat- - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VarValue SclPrd = 0;// scalar product
  for( Index i = 0 ; i < NumVar ; i++ )
   SclPrd += dir[ i ] * ( LambdaHat[ i ] - Lambda[ i ] );

  HatEpsilon += alpha * ( FiLambda - FiHat + SclPrd );

  // compute the linearization error SigmaHat of Gi[wFi] in LambdaHat  - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SclPrd = 0;
  for( Index i = 0 ; i < NumVar ; i++ )
   SclPrd += Gi[ i ] * ( Lambda[ i ] - LambdaHat[ i ] );

  SigmaHat += (1.0 - alpha ) * ( FiHat - FiLambda + SclPrd );

  // compute the LambdaHat and FiHat - - - - - - - - - -  - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for(Index i=0; i<NumVar;i++){
	  VarValue g=LambdaHat[i];//lambdahat_{i-1}
	  LambdaHat[i]=Lambda[i]*alpha + (1.0- alpha)*g;
  }

  FiHat = alpha * FiLambda + (1.0 - alpha) * FiHat;// f(lambdahat)=alpha*f(lambda)+(1-alpha)*f(lambdahat)
  //f is a linear C05Function

  // compute the linearization error HatEpsilon of Dir in LambdaHat  - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  HatEpsilon = alpha * SigmaHat + ( VarValue(1) - alpha ) * HatEpsilon;
  }  // end( alpha < 1 )- - - - - - - - - - - - - - - - - - - - - - - - - - -

 //VLog sul valore di FiHat

 }  // end( SubGradSolver::ChgHatlambda() )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- End Class SubGradSolver------------------------*/
/*--------------------------------------------------------------------------*/

/*-------------------------- CLASS Stepsize --------------------------------*/
/*--------------------------------------------------------------------------*/
/** These methods provide derived classes of Stepsize (who are not friends of
 *  SubGrad, while Stepsize is) access to protected data of SubGrad.
 *--------------------------------------------------------------------------*/

VarValue Stepsize::GetCoeffDefl( void ){
 return( Solver->alpha );
 }

/*--------------------------------------------------------------------------*/

VarValue Stepsize::GetGiNorm( void ){
 return( Solver->NrmGi );
 }

/*--------------------------------------------------------------------------*/

VarValue Stepsize::GetDNorm( void )
{
 return( Solver->NrmDir );
 }

/*--------------------------------------------------------------------------*/

VarValue Stepsize::GetdGk( void )
{
 return( Solver->dGk );
 }

/*--------------------------------------------------------------------------*/

VarValue Stepsize::GetdkM1Gk( void )
{
 return( Solver->dM1Gk );
 }

/*--------------------------------------------------------------------------*/

Index Stepsize::GetNItIcr( void )
{
 return( Solver->NItIncr );
 }

/*--------------------------------------------------------------------------*/

VarValue Stepsize::ReadFkVal( void )
{
 return( Solver->FiLambda );
 }

/*--------------------------------------------------------------------------*/

VarValue Stepsize::ReadFiBar( void )
{
 return( Solver->FiBar );
 }

/*--------------------------------------------------------------------------*/

VarValue Stepsize::GetGlobalLipschitz( void ){
 	return(Solver->f_lf->get_Lipschitz_constant());
 }

/*--------------------------------------------------------------------------*/

VarValue Stepsize::ReadBestFiVal(void){
  return(Solver->FiBest);
}

/*--------------------------------------------------------------------------*/

VarValue Stepsize::GetLowerbound(void){
  return(Solver->Lowerbound);
}

/*--------------------------------------------------------------------------*/

VarValue Stepsize::GettStar(void){
  return(Solver->tStar);
}

/*--------------------------------------------------------------------------*/

void Stepsize::set_par(const idx_type par , const int value ){
		throw( std::invalid_argument( "It is not a valid parameter" ) );
}

/*--------------------------------------------------------------------------*/

void Stepsize::set_par(const idx_type par , const double value){
		switch(par){
			case(dblLpsFct):
				LpsFct=value;
			break;
		default:
			throw(std::invalid_argument("Unknown parameter"));
		}
}

/*--------------------------------------------------------------------------*/

int Stepsize::get_int_par(const idx_type par) {
		throw( std::invalid_argument( "There is not any int parameter" ) );
}

/*--------------------------------------------------------------------------*/

double Stepsize::get_dbl_par(const idx_type par) {
		switch(par){
			case(dblLpsFct):
				return(LpsFct);
				break;
			default: throw(std::invalid_argument("Unknown parameter"));
		}	
}

/*--------------------------------------------------------------------------*/
	
int Stepsize::get_dflt_int_par( const idx_type par ) {
		throw( std::invalid_argument( "There is not any int parameter" ) );
 }

 /*--------------------------------------------------------------------------*/
 
 double Stepsize::get_dflt_dbl_par( const idx_type par ) {
		switch(par){
			case(dblLpsFct):
				return(dflt_dbl_par[dblLpsFct]);
			break;
		default:
			throw(std::invalid_argument("Unknown parameter"));
    }
}

/*--------------------------------------------------------------------------*/

Stepsize::StepsizeFactoryMap & Stepsize::f_factory( void ){
	static StepsizeFactoryMap s_factory;
	return( s_factory );
	}

/*--------------------------------------------------------------------------*/
/*---------------------------- end Stepsize.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Deflection ------------------------------*/
/*--------------------------------------------------------------------------*/
/** These methods provide derived classes of Deflection (who are not friends
 *  of SubGradSolver, while Deflection is) access to protected data of 
 *  SubGradSolver.
 *--------------------------------------------------------------------------*/

VarValue Deflection::GetStepsize( void ){
 return( Solver->step );
 }

/*--------------------------------------------------------------------------*/

VarValue Deflection::GetGiNorm( void ){
 return( Solver->NrmGi );
 }

/*--------------------------------------------------------------------------*/

VarValue Deflection::GetDNorm( void ){
 return( Solver->NrmDir );
 }

/*--------------------------------------------------------------------------*/

VarValue Deflection::GetdGk( void ){
 return( Solver->dGk );
 }

/*--------------------------------------------------------------------------*/

VarValue Deflection::GetSigma( void ){
 return( Solver->Sigma );
 }

/*--------------------------------------------------------------------------*/

VarValue Deflection::GetEpsilon( void ){
 return( Solver->Epsilon );
 }

/*--------------------------------------------------------------------------*/

VarValue Deflection::ReadFVal( void ){
 return( Solver->FiLambda );
 }

 /*--------------------------------------------------------------------------*/

VarValue Deflection::GetGlobalLipschitz( void ){
 	return(Solver->f_lf->get_Lipschitz_constant());
 }

/*--------------------------------------------------------------------------*/

VarValue Deflection::ReadFiVal( void ){
    return(Solver->FiBar);
 }

/*--------------------------------------------------------------------------*/

Index Deflection::GetNrIter(void){
  return(Solver->ParIter);
}

/*--------------------------------------------------------------------------*/

VarValue Deflection::ReadBestFiVal( void ){
 //if( wFi!= Inf<Index>() )
  //throw( std::invalid_argument( "ReadBestFiVal: wFi != INF" ) );
	//if(IsConvex)
 		return( Solver->FiBest );
 	//else {
 	//	if(IsConcave)
 	//		return (-FiBest);
 	//	else
 	//		throw (std::invalid_argument("ReadBestFiVal: Fi() is neigther convave or convex"));
	}

/*--------------------------------------------------------------------------*/

double Deflection::GetEpsLin(void){
  return(Solver->EpsLin);
}

/*--------------------------------------------------------------------------*/

Deflection::DeflectionFactoryMap & Deflection::f_factory( void ) {
static DeflectionFactoryMap s_factory;
return( s_factory );
}

/*--------------------------------------------------------------------------*/
/*------------------------- end(Deflection.cpp) ----------------------------*/
/*--------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/*--------------------------End File SubGradSolver.cpp-----------------------*/
/*---------------------------------------------------------------------------*/
