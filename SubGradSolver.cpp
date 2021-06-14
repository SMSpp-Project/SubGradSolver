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
#include "PrimalDualClass.h"

#include "PrimalDual.h"
#include "FumeroTV.h"
#include "ColorTV.h"
#include "Polyak.h"
#include "Volume.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/
#define LOG_LEVEL 0
#if (LOG_LEVEL>1)
  #define VLOG( l , x ) if( LogVerb > l ) *f_log << x 
  #define VLOG2( l , c , x ) if( ( LogVerb > l ) && c ) *f_log << x
#else
  #define VLOG( l , x )
  #define VLOG2( l , c , x )
#endif

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
			{ "intSGPar1" , "intSGPar3" , "intSGPar4","intSGPar5",
			 "intRstAlgPar" };

// define and initialize here the vector of double parameters names
const std::vector< std::string > SubGradSolver::dbl_pars_str =
			{ "dblSGPar2" ,"dbltStar", "dblEpsLin","dblEInit",
			"dblEFnalt","dblEDcrs","dblEStps" };	
			
// define and initialize here the vector of string parameters names
const std::vector< std::string > SubGradSolver::str_pars_str ={
"Stepsize_name_par", "Deflection_name_par"};		

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
			 
// define and initialize here the map for string parameters names
const std::map< std::string , SubGradSolver::idx_type > SubGradSolver::str_pars_map =
		  { {"Stepsize_name_par", SubGradSolver::Stepsize_name_par},
		    {"Deflection_name_par", SubGradSolver::Deflection_name_par}	};

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

// define and initialize here the default string parameters
const std::vector<int> SubGradSolver::dflt_str_par ={
 "Volume",
 "null"
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

 if( f_Block == block )  // registering to the same Block
  return;                // cowardly and silently return

 if( f_Block ) {  // changing from a previous oracle - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  guts_of_destructor();   // deallocate memory
  }
 Solver::set_Block( block );  // attach to the new Block

 if( ! f_Block )  // that was actually clearing the Block
  return;         // all done

  bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
                       "SubGradSolver: unable to lock the Block" ) );

 /* Two types of block can be handled by the SubGradSolver:

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

 if( sb.empty() ) { //no sub-Block
  // the objective function of the block must be a C05Function  - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
  if( ! obj )
   throw( std::logic_error( "objective is not a FRealObjective" ) );

  auto c05f = dynamic_cast< C05Function * >( (obj)->get_function() );
  if( ! c05f )
   throw( std::logic_error( "the objective is not a C05Function" ) );

  IsConvex = c05f->is_convex();
  if( ( ! IsConvex ) && ( ! c05f->is_concave() ) )
   throw( std::invalid_argument(
			     "only convex or concave objectives allowed " ) );

  if( ( IsConvex && ( obj->get_sense() == Objective::eMax ) ) ||
      ( ( ! IsConvex ) && ( obj->get_sense() == Objective::eMin ) ) )
   throw( std::invalid_argument( "can only minimize convex / maximize concave"
				 ) );  

  v_c05f.push_back( c05f );
  f_lf = nullptr;
  }
 else { //there are sub-Block
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
  //printf ("nrfi=%zu\n", sb.size());
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
  // continue;
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
  
  // if some Variable are present, they are of the ColVariable type - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  
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

 /* auto v_d_Variable = f_Block->get_dynamic_variables();
 for( auto & el : v_d_Variable ) {
  if( un_any_thing_0( std::list< ColVariable > , el , ++NumBVar ) )
   continue;
  if( un_any_thing_1( std::list< ColVariable > , el ,
		      NumBVar += var.size() ) )
   continue;
  if( un_any_thing_K( std::list< ColVariable > , el ,
		      NumBVar += var.num_elements() ) )
   continue;
  throw( std::logic_error( "some dynamic Variable is not a ColVariable" ) );
  }
  */
 if( NumBVar < NumVar )
  throw( std::logic_error( "too few ColVariable in the Block" ) );

 // check that the Variable in the Block agree with that in the C05Function- -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::vector<ColVariable *> LamBVcblr( NumBVar );

 Index cnt = 0;
 for( auto & el : v_s_Variable )
  un_any_static( el , [ & ]( ColVariable & sv ) { LamBVcblr[ cnt++ ] = & sv;
                  } , un_any_type<ColVariable>() );

  /* for( auto & el : v_d_Variable )
  un_any_dynamic( el , [ & ]( ColVariable & sv ) { LamBVcblr[ cnt++ ] = & sv;
                  } , un_any_type<ColVariable>() );
*/
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
 //printf("NrFi=%d\n", NrFi);
 
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
 
 /* if( ( SGPar1 & PrjDIRM1 ) && ( SGPar1 & PrjSG ) )
   SGPar1 &= ~PrjDIR;*/
   
  // reset algorithm  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
  // checking for 0-th component - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 //Index * SG0Bse;
  //Index SG0BDm =GetGi( gi , SG0Bse ,0, NrFi);
  //printf("SG0BDm=%d\n", SG0BDm);
    for(Index i=0; i<NumVar; i++){
      //Gi[i]=gi[i];
      Gi[i]=0;
      //printf("Iter n %d\n", i);
    }
  //ZeroComp = SG0BDm > 0;
  ZeroComp=true;
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
	//if(stepsize)
	InitStepsize();
	
	InitDeflection();  

  // other initializations- - - - - - - - - - - - - - - - - - - - - - - - - -
  fs = ThinComputeInterface::kUnEval;
 
  EmptySet = false;       // the feasible set is not empty

/*-----------------------------------------------------------------------------*/
	/*for(Index i=0; i<NrFi; i++){
		IsConvex = IsConvex && (v_c05f[i]->is_convex()); //is true iff every component of the objective function is convex
		IsConcave = IsConcave && (v_c05f[i]->is_concave()); //is true iff every component of the objective function is concave
	}*/
	
	UpdtLowerBound();
	
	AllocateUBLB(0, NumVar);

 //ReSetAlg( RstAlgPar );  // Lambda is reset inside
 if( ! owned )
  f_Block->unlock( f_id );

 //printf("Finish set_block\n");
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
 Result = Solver::kUnEval; //risultato inizializzato a default
 SCalls++; //number of calls to compute() included this one
 Index comp;	//a component
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 do {
  // function log: print the information relative to the stability center
  // and the best point - - - - - - - - - - - - - - - - - - - - - - - - - - -
printf("\n\nCiclo nuovo\n\n");
  Log1();
 //printf("\n Ciclo %d\n", ParIter);
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
  if(FiBar < Inf<VarValue>()){ 
    printf("Sono nell'if, FiBar<inf\n");
	//if the function value (at Lambdabar) is unknown  skip the direction and stepsize
   // computation and go straight to FiAndGi()
   //printf("Ciclo interno\n");
   // direction log: print the information relative to DIRM1  - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Log2();

   // FormD call: a new direction is provided, but Dir is not updated
   // yet because AddVariables() and RemoveVariables() need to work
   // with the prior vector Dir - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   FormD( );
  printf("Sono fuori da formD()\n");
	// check the status and take the action, only during
   // outer iterations  - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //printf("Log4\n");
  printf("Inniter after FormD()=%d\n", InnIter);
   if( !InnIter ){
    for(Index i=0; i<NrFi; i++){
      fs=v_c05f[i]->compute(false);
    }
    printf("Ricompute\n");
   }

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
   
   //: note that the previous DIR is not saved
   // as DirM1, but it is merely discharged; only and one vector is kept
   // for storing the direction- - - - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
   SaveDir();
   printf("Sono fuori da savedir\n");
   //printf("Log5\n");
   // if the Stepsize object needs the scalar product between DirM1 and Gi (<d_{i-1}, g_i>),
   // then the direction vector is kept- - - - - - - - - - - - - - - - - - - -

   if( stepsize->NeedsdkM1Gk() && ( ! deflection ) ){
     printf("Needsdkm1gk=%d\n", stepsize->NeedsdkM1Gk());
	for(Index i=0; i<NumVar; i++)
		dir[i]=Gi[i];
   }
   // direction log: print the information relative to DIR   - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for(Index i=0; i<NumVar;i++)
    printf("After formd e savedir dir[%d]=%f\n", i, dir[i]);
   Log2();
	
	if( step <= 1e-8 * tStar ) { // it'a a shortstep
    printf("SmallStep\n");
    CSmallStep++; // another shortstep 
    if( CSmallStep > ( ( NItIncr > 0 )? 100 * NrFi : 100 ) ) { // If the number of incremental iteration after an outer iteration (NItIncr) is >0 then I compare CSmallStep with 100*NrFi (numbers of component of the obiective sum function), else with 100
	// If the number of consecutive shortstep is higher than the limit value, then I stop the iteration and return the status
     Result = kStopIter;
     break;
     }
    }
   else{// it's a normal step 
   printf("normal step=%f\n", step);
    CSmallStep = 0;//The counter is reinizialized because we only check about consecutive short steps
   }  // end( search direction computation, when there are no info about the function valuated in LambdaBar) - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  }




  // check for optimality  - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //printf("\n Inniter=%d IsOptimal=%d\n", InnIter, IsOptimal());
  if( ( ! InnIter ) && IsOptimal() ){
    printf("\n entro qui\n");
	   // run optimality events
   int res = eContinue;
   for( auto & ev : v_events[ eBeforeTermination ] )
    if( ( res = ev() ) != eContinue )
     break;

   if( res == eStopError ) {
    VLOG( 1 , " ~ stop (error) by optimality event" << std::endl );
    Result = kError;
    break;
    }
    
	if(res==eForceContinue){
		throw(std::logic_error("This is not allowed"));
		break;
	}
   VLOG( 1 , " ~ stop (optimal)" << std::endl );
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
   VLOG( 1 , " ~ stop due to time-periodic event" << std::endl );
   break;
   }
  // calculate the tentative point Lambda1 (not projected yet)   - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  FormLambda1();
  VLOG( 1 , std::endl << "           " << " |Lambda1| = " << Norm( Lambda , NumVar )  );

  // a real iteration (iterations where Fi() is not evaluated count as well) -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ParIter++; //new iteration
  printf("pariter=%d\n", ParIter);
  // is an inner iteration?- - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  printf("NIrIncr=%d GiEvaltns=%d\n", NItIncr, GiEvaltns);
  InnIter = NItIncr? GiEvaltns % ( NItIncr + 1 ): 0; 
  //InnIter takes GiEvaltns%(NItIncr+1) if NItIncr!=0, else 0
  //GiEvaltns is the number of calls to Gi()
  printf("InnIter calcolato=%d\n", InnIter);

  LHasProj = false; //lambda has not been projected yet in this iteration
  
  // pick up component comp (defined at the beginning of the function) - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//printf("Inniter=%d\n", InnIter);
  if (InnIter){
	//	while (!comp)// I don't want the first component
			comp=myrandom(NrFi);//I take one of the component
	}
	else{
    printf("Inniter=0\n");
		comp=Inf<Index>();
  }
    printf("comp=%u\n", comp);
	// project Lambda1 onto the feasible set - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
	while (1){
    //printf("\n In while(1)\n");
    bool t=ProjFsb(Lambda, EmptySet);
   printf("Risultato ProjFsb=%d\n", t);
		if(t){//check if lambda components were changed in the projection
			LHasProj=true;
      LHasChgd=true;
			VLOG( 1 , "\n ~ |Lambda1| = " << Norm( Lambda , NumVar )  );
			if(EmptySet)
				break;
			//continue while the set is feasible: so it will be possible to take a point	
			}

		// never evaluate the objective function at the same point - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		if( ! LHasChgd )//if the point lambda is the same as the previous iteration
			throw (std::logic_error("SubGradSolver::compute() this should not happen"));
	
		// compute Fi( Lambda1 ) and a subgradient / constraint Gi[wFi]- - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
		FiAndGi( comp );//update LmbdBst
    for(Index i=0; i<NumVar; i++)
    printf("Gi[%d]=%f\n", i ,Gi[i]);
		// if the Lambda1 is infeasible, go back to the projection - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   printf("FuoridaFiandGi\n");
   printf("FiLambda fuori da fiandgi=%f\n", FiLambda);
		if( FiLambda < Inf<VarValue>() )
			break;
    
   			} // end function evaluation - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// project the subgradient Gi[wFi] over the tangent cone at tentative point
   // Lambda1- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( SGPar1 & PrjSG ){// check if Gi[wFi] is not projected yet
  printf("SGPAR1 & PrjSG projt\n");
   ProjTng( Gi , Lambda );
  }
  NrmGi = 0;//Norm( Gi , NumVar );//compute the sqared norm of Gi[wFi]
  for(Index i=0; i<NumVar; i++)
  NrmGi=NrmGi + (Gi[i]*Gi[i]);
  printf("NormGi afer eventually projection=%f\n", NrmGi);
  // compute the scalar product between Dir and Gi[wFi] at Lambda1- - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -	
  if( stepsize->NeedsdkM1Gk() ) {
  printf("needsdmkm1gk\n");
   dM1Gk = ScalarProduct( dir , Gi , NumVar );
   dM1GkDone = true;
   }
  else
   dM1GkDone = false;
    // print the logs regarding Fi and Gi[wFi] at Lambda1 - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   //log su fi(lambda1) e gi(wfi) 
   printf("Scalar product=%f\n", ScalarProduct( dir , Gi , NumVar ));
   // check for feasibility- - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
   if( EmptySet ) {
     //printf("kInfeasible\n");
		Result = kInfeasible;
   break;
   }
   
   if (FiLambda==-Inf<VarValue>()){
      //printf("kUnbounded\n");
	   Result=kUnbounded;
	   break;
   }
  // check whether or not the Lower Bound has changed - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //printf("before\n");
	UpdtLowerBound();
printf("Updatelb\n");
  VLOG2( 1 , Lowerbound > - Inf<VarValue>() , "\n ~ UB = " << Lowerbound );
   
   // if FiBar (the function at LambdaBar) is unknown, or no deflection is
  // used, then move to Lambda1 - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  printf("FiBar=%f\n", FiBar);
   if( ( FiBar == Inf<VarValue>() ) || ( ! deflection ) ) {
		GotoLambda1();
    printf("Go to lambda 1\n");
		DoSS = true;//A serious step is done
		VLOG( 1 , std::endl );
		continue;
   }
   printf("non eseguo go to lambda 1\n");
   // Null Step (NS) vs Serious Step (SS)- - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  DoSS = deflection->DoSS();//perchéè sempre falso?
 printf("In SubGrad:doSS=%d\n", deflection->DoSS());
  if( DoSS ) {  // SS- - - - - - - - - - - - - - - -
  printf("In SubGrad:delta ss=%f\n", deflection->Delta());
   VLOG2( 1 , deflection->Delta() < Inf<VarValue>() , std::endl << " SS[" << CSSCntr );

   GotoLambda1();
  
  CSSCntr++;//new serious step is done
  CNSCntr = 0;// turn the number of consecutive null steps to 0
  printf("CSSCntr=%d\n", CSSCntr);
  }else {  // NS- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   VLOG2( 1 , deflection->Delta() > - Inf<VarValue>() , std::endl << " NS[" << CNSCntr );
  printf("In SubGrad:delta ns=%f\n", deflection->Delta());
   UpdateSigma();//need to modify the error 
  printf("Sigma=%f\n", Sigma);
   CNSCntr++;
   CSSCntr = 0;
   printf("CNSCntr=%d\n", CNSCntr);
   }  // end SS/NS - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  VLOG
  ( 1 , "]: DFi (" << FiBar - FiLambda << ") < m1 * Dv ("
	     << deflection->Delta() << ")" );
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
   printf("nuovociclo\n");
  } while( ( ! MaxIter ) || ( ParIter < MaxIter ) );
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //printf("Return=%d\n", Result);
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
		SGPar1 = value;
    if( ( SGPar1 & PrjDIRM1 ) && ( SGPar1 & PrjSG ) )
      SGPar1 &= ~PrjDIR;
    if( ! deflection )
      SGPar1 &= ~( PrjDIRM1 | PrjDIR );
	break;
   case(intSGPar3):
		SGPar3 = value;
    if( ! deflection )
      SGPar3 |= StepRST;
	break;
	case(intSGPar4):
    if(value==0 || value==1){
		  if( NItIncr ) // forbidden for incremental version
        SGPar4 = 0;
      else
        SGPar4 = value;
    }
    else 
      throw(std::invalid_argument("SGPar4 value not allowed"));
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
				SGPar2 = value;
        InnIter=false;
        ZeroComp=true;
        if (ZeroComp)
           NItIncr = ceil( double( NrFi + 1 ) * SGPar2 );
        else
           NItIncr = ceil( double( NrFi ) * SGPar2 );
        printf("In set_par NItInctr=%d\n", NItIncr);
        if( NItIncr ) {
           SGPar4 = 0;
           if( deflection )
              throw( std::invalid_argument( "set_par: deflection is not allowed" ) );
        }
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
			/*case(dblEpsLin):
				EpsLin=value;
			break;*/
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

void SubGradSolver::set_par( idx_type par , std::string && value  ){
	switch(par){
		case(Stepsize_name_par):
			if(value=="FumeroTV" || value=="ColorTV" || value=="Polyak")
				Stepsize_name=value;
			else
				throw(std::invalid_argument("Stepsize method not allowed"));
			break;
		case(Deflection_name_par):
			if(value=="Volume" || value=="PrimalDual" || value=="null")
				Deflection_name=value;
			else
				throw(std::invalid_argument("Deflection method not allowed"));
			break;
		default: CDASolver::set_par(par, value);
		}
} //end (SubGradSolver::set_par(string))




/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void SubGradSolver::set_log( std::ostream *log_stream){
 f_log=log_stream;
  if( LogVerb > 1 ) {
   *f_log << std::endl << "SubGradSolver: Vars = " << NumVar << "NrFi="<< NrFi <<" ("
	   << ") ~  NItIncr = " << NItIncr <<  " ~ EpsLin = "
	   << EpsLin << " ~ Einit = " << EInit << " ~ t* = " << tStar;
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
		/*case(dblEpsLin):
			return(EpsLin);
		break;*/
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
printf("InIsOptimal:\n  tStar=%f, EpsLin=%f, Epsilon=%f, NrmDir=%f, FiBest=%f\n", tStar, EpsLin, Epsilon, std::sqrt( NrmDir ), FiBest);
 if( FiL == Inf<VarValue>() )
  return( false );
 else {
  if( FiL < 0 ) FiL = - FiL;
  if( FiL < 1 ) FiL = 1;  //max{1, |f_i^{rec}|}

  if( eps <= 0 )
   eps = EpsLin;
  double d=0.1;
  if( SGPar4 )
   return( tStar * std::sqrt( NrmDir ) + std::min( Epsilon , HatEpsilon )
	   <= eps * FiL );
  else
  return( tStar * std::sqrt( NrmDir ) + d <= eps * FiL );
  // return( tStar * std::sqrt( NrmDir ) + Epsilon <= eps * FiL );
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
	//FiBest=FiBar=Inf<VarValue>();
	
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
if ((SGPar1 & PrjSG) && (! DoSS) ) //status check(fs==//status check)) 
//g_i is already projected but no serious step has been done
//The projection has to be recomputed
{
  printf("DoSS in FormD =%d\n", DoSS);
	ProjTng(Gi, LambdaBar);
NrmGi=Norm( Gi, NumVar);// the squared norm
printf("Normgi in FormD=%f\n", NrmGi);
dM1GkDone = false;//the scalar product \f$d_{i-1}^{\top} g_i\f$ has to be computed
  }
// subgradient log: print the information relative to Gi[wFi] - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Log4();
//printf("F1\n");
 // project DirM1 and compute its norm - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (SGPar1 & PrjDIRM1) //status check(fs==//status check))
	//if d_{i-1} has been projected
	{if(SGPar1 & PrjDIR){
		//if d_i has been projected over the tangent cone of the previous Lambdabar
		if (DoSS){ //if a serious step has been done
			ProjTng(dir, LambdaBar);// mandatory to reproject
			NrmDir=Norm(dir, NumVar);
          printf("In FormD1 DoSS=%d NrmDir=%f\n", deflection->DoSS(), NrmDir);
			dM1GkDone = false; //the scalar product \f$d_{i-1}^{\top} g_i\f$ has to be computed
		}	
	}
	else{//is needed to project d_i
    printf("Not projected yet\n");
		ProjTng(dir, LambdaBar);
		NrmDir=Norm(dir, NumVar);
    printf("In FormD NrmDir=%f\n", NrmDir);
    for (Index i=0; i<NumVar; i++)
        printf("dir[%d]=%f\n", i, dir[i]);
		dM1GkDone = false; //the scalar product \f$d_{i-1}^{\top} g_i\f$ has to be computed
	}
}	

// (re)compute the scalar product between DirM1 and Gi[wFi] at Lambda - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (deflection){
	if(! dM1GkDone)//if the scalar product has not been already computed
		dM1Gk=ScalarProduct( Gi, dir, NumVar);
  printf("deflection!=null dM1Gk=%f\n", dM1Gk);
}
else{
	//no deflecton is done -> d_i=g_i
	dGk=NrmGi;
  NrmDir=NrmGi;// the scalar product <d_i, g_i>
	if(! dM1GkDone)
    dM1Gk=dGk;// <d_{i-1}, g_i>=<d_i, g_i>
}

dM1GkDone=true; //the scalar product <d_{i-1}, g_i> has been done 
//printf("F3\n");
/*-------------------------------------------------------------------------------*/
/* Check on the approach: stepsize-restricted or deflection-restricted.
This information is coded bit wise in SGPar3: second bit is equal to 1
in the first approach, equal to 0 in the secon one*/
if (SGPar3 & StepRST){
  printf("Stepsize restricted\n");
	//stepsize-restricted approach: deflection, compute, stepsize
	if (deflection){
		//deflection computation
    printf("deflection computation\n");
		deflection->NewDEF();// inizialization
		alpha=deflection->GetDFLCoeff();// get alpha 
    printf("alpha=%f\n", alpha);
	}
VLOG( 1 , " ~ alpha = " << alpha );
}
else{
  printf("Deflection restricted\n");
	//deflection-restricted approach: stepsize, deflection, compute
	dGk=dM1Gk;// is not important the distinction between them
	printf("InnIter in deflection restricted =%d\n", InnIter);
	if( ! InnIter){// if is an outer iteration is needed to recompute stepsize
	stepsize->NewStep();
  printf("NewStep\n");
  }
	
	step=stepsize->GetStepsize(InnIter);// get the stepsize according to the type of iteration
  printf("Step dopo getstepsize =%f\n", step);
	
 VLOG( 1 , " ~ step = " << step );
	
	deflection->NewDEF();// the case of NO DEFLECTION is not allowed
	alpha=deflection->GetDFLCoeff();// get alpha
 printf("In deflection restricted alpha=%f\n", alpha);
	
 VLOG( 1 , " ~ alpha = " << alpha );
	/*-------------------------------------------------------------------------------------------*/
	
	/* Check if any safe rule is taken according to SGPar3 coded bit-wise: 
	the first bit is equal to 1 if the rule is used, it's equal to 0 otherwise*/
	
	if((SGPar3 & SafeRULE) && !InnIter){
    printf("Saferule\n");
		//if a safe rule is used and there is an outer iteration
		VarValue target= (step * NrmDir)/(FiBar-stepsize->GetLev()+ (step + NrmDir)); //(v_i*||d_{i-1}||^2)/(f_i-f_i^lev + v_i*||d_{i-1}||^2 )
		//f_i^lev=f_i^rec-\delta_i
		//f_i^rec==min{f_k:k=1,..,1}

		//\delta_i>0 depends on the choice of vanish or nonvanishing method
		alpha=std::min(1.0, std::max(alpha, target));// to respect the right inequality for the safe rule
    printf("Safe rule target=%f e alpha=%f\n", target, alpha);
		
   VLOG( 1 , " ~ fixed alpha = " << alpha );
	}
}
for (Index i=0; i<NumVar; i++)
printf("dir[%d]=%f in formd\n", i, dir[i]);

// print Dir norm - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
VLOG( 1 , " ~ trial |D| = " << NrmDir );

}// end SubGradSolver::FormD()

/*--------------------------------------------------------------------------*/
void SubGradSolver::SaveDir(void){
// compute the linearization error of Dir at LambdaBar- - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //printf("alpha in Savedir=%f\n", alpha);
 if(alpha==1){
  printf("savedir alpha=1\n");
  Epsilon=Sigma;
 }

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
	 printf("Savedir deflection\n");
	 if(alpha==1){
		 for(Index i=0; i<NumVar; i++){
			dir[i]=Gi[i];
      printf("dir[%d]=%f\n", i, dir[i]);
     }
	 }
	 else{
		 for(Index i=0; i<NumVar; i++){
       printf("Savedir alpha !=1\n");
			 VarValue g=dir[i];//if we do the k-th iteration this is d_{k-1}[i]
			 dir[i]=alpha*Gi[i]+(VarValue(1)-alpha)*g;
       //printf("dir[%d]=%f\n", i, dir[i]);
		 }
			
	 }
	 for (Index i=0; i<NumVar; i++)
  printf("dir[%d]=%f in savedir1\n", i, dir[i]);
	 
	 // project Dir on tangent cone at LambdaBar  - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if(SGPar1 & PrjDIR){// capire perchè fa questo controllo
		ProjTng(dir, LambdaBar);
    printf("Projdir in savedir\n");
  }
	
	// compute the (squared) norm of the direction and its scalar product with
  // the subgradient   - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	NrmDir = 0;//Norm( Gi , NumVar );//compute the sqared norm of Gi[wFi]
  for(Index i=0; i<NumVar; i++)
  NrmDir=NrmDir + (dir[i]*dir[i]);
  printf("NrmDir= in safedir %f\n", NrmDir);
	dGk = ScalarProduct( Gi , dir , NumVar );
  printf("dgk in savedir=%f\n", ScalarProduct( Gi , dir , NumVar ));

  }  // end( deflection )  - - - - - - - - - - - - - - - - - - - - - - - - - -
	
	// stepsize computation - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( SGPar3 & StepRST ) {//check on the approach: stepsize restricted
  if( ! InnIter ) // newstep only in outer iterations
   stepsize->NewStep( );

  step = stepsize->GetStepsize( InnIter );
  
   VLOG( 1 , " ~ step = " << step );
  // safe rule - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( SGPar3 & SafeRULE ) {//check on the safe rule
   stepsize->SetMaxBeta( alpha );//capire
   step = stepsize->GetStepsize( InnIter );
    VLOG( 1 , " ~ (fixed) step = " << step );
   
   }

  }  // end stepsize-restricted scheme - - - - - - - - - - - - - - - - - - - -
  for (Index i=0; i<NumVar; i++)
printf("dir[%d]=%f fine di savedir\n", i, dir[i]);
 }  // end( SubGradSolver ::SaveDir() )

/*--------------------------------------------------------------------------*/
void SubGradSolver::FormLambda1(void){
	if(deflection){
     for (Index i=0; i<NumVar; i++)
  printf("dir[%d]=%f in FormLambda1\n", i, dir[i]);
    printf("\n deflection: step= %f\n", step); 
		for( Index i = 0 ; i < NumVar ; i++ ){
			Lambda[i]=LambdaBar[i] - step * dir[i];
      LamVcblr[ i ]->set_value( Lambda[ i ] );
      //printf("e dir=%f\n", dir[i]);
    }
	}
	else{
     for (Index i=0; i<NumVar; i++)
    printf("Gi[%d]=%f FormLambda1\n", i, Gi[i]);
    printf ("\n no deflection: step= %f", step);
		for( Index i = 0 ; i < NumVar ; i++ ){
			Lambda[ i ] = LambdaBar[ i ] - step * Gi[ i ];
      LamVcblr[ i ]->set_value( Lambda[ i ] );
      //printf (" e dir=Gi=%f\n", Gi[i]);
    }
	}
 	LHasChgd = true;// lambda has changed since the last call to FiandGi()
   for (Index i=0; i<NumVar; i++)
    printf("FormLambda1 LamVcblr[%d]=%f Lambda[%d]=%f\n",i, LamVcblr[i], i , Lambda[i]);
}	//end(SubGradSolver::FormLambda1())

/*--------------------------------------------------------------------------*/
void SubGradSolver::FiAndGi(  Index wFi ){
 // pass Lambda1 - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( LHasChgd ) { 
   printf("In FiandGi Lambdahaschngd\n"); // if Lambda has changed, set the new one 
  SetLambda( Lambda );
  LHasChgd = false;
//for(Index i=0; i<NumVar; i++)
//printf("Lambdabar[%d]=%f\n", i, LambdaBar[i]);
  }
  for(Index i=0; i<NrFi; i++){
    fs=v_c05f[i]->compute(false);
    if(fs!=ThinComputeInterface::kOK)
      throw(std::logic_error("Function compute fails"));
  }
   FiEvaltns++;
   //printf("After compute fs=%d\n", fs);
 // if outer iteration, evaluate Fi in Lambda1  - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( wFi == Inf<Index>() ) { 
   printf("Caso comp=inf in fiandgi\n");
  FiLambda=v_c05f[0]->get_value(); // call Fi()- - - - - - - - - - - - - - - - - -
  for(Index i=1; i<NrFi; i++)
    FiLambda = FiLambda+ v_c05f[i]->get_value();//f_lf->get_value();//f(lambda)
    //printf("FiLambda=%f\n", FiLambda);
 
  if( FiLambda == -Inf<VarValue>() )  // Fi() unbounded (-inf)
   return;
  // update FiBest, if necessary - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( FiLambda < FiBest ) {
   FiBest = FiLambda;
   KpBstL=true;
    for(Index i=0; i<NumVar; i++)
		LmbdBst[i]=Lambda[i];
   }
  }  // end( function evaluation) - - - - - - - - - - - - - - - - - - - - - - -

 // if Lambda is unfeasible, a new subgradient is computed; otherwise a new
 // constraint is added - - - - - - - - - - - - - - - - - - - - - - - - - - -

if( FiLambda == Inf<VarValue>() ){// it is a constraint
	throw(std::logic_error("It's out of the feasible region"));
}
else{
  EvalGi( wFi );  // the item is a subgradient 
 }
 
}//end(SubGradSolver::FiandGi())

/*--------------------------------------------------------------------------*/
void SubGradSolver::GotoLambda1( void ){
 // compute the linearization error of Dir at the new LambdaBar  - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
printf("Inniter in gotolamda1=%d\n", InnIter);
 if( FiBar < Inf<VarValue>() && deflection ) {
   printf("Gotolambda1\n");
  VarValue SclPrd;
  if( ! LHasProj ){
    printf("LHasProj=true\n");
   SclPrd = step * NrmDir;
  }
  else {
    printf("Lambdahasnotproected\n");
   SclPrd = 0;
   for( Index i = 0 ; i < NumVar ; i++ )
    SclPrd += dir[ i ] * ( LambdaBar[ i ] - Lambda[ i ] );
   }
  printf("Scalarprod=%f\n", SclPrd);
  Epsilon += FiLambda - FiBar + SclPrd;
  printf("Epsilon=%f\n", Epsilon);
  }

 // move Lambda to LambdaBar= - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	for(Index i=0; i <NumVar; i++)
		LambdaBar[i]=Lambda[i];
	if( ! InnIter ){
	FiBar = FiLambda;
  printf("Goto FiBar=%f\n", FiBar);
  }
	
  for(Index i=0; i<NumVar; i++)
  printf("GotoLambda1: LambdaBar[%d]=%f\n", i, LambdaBar[i]);
 }  // end( SubGradSolver::GotoLambda1() )

/*--------------------------------------------------------------------------*/
Index SubGradSolver::GetGi(VarValue *sub, Index * sgb, Index strt, Index stp){
	Index NumFun = v_c05f[ 0 ]->get_num_active_var();
  //printf("Sono dentro Gi\n");
	if( stp > NumFun )
		stp = NumFun;
		auto range=std::make_pair(strt, stp);
//	HasLin= v_c05f[strt]->has_linearization(true);
  //HasLin=true;
  for(Index k=strt; k<stp; k++){
  HasLin= v_c05f[k]->has_linearization(false);
// HasLin=true;
  //printf("HasLin=%d\n", HasLin);
	if(HasLin){
    //printf("Linearizz\n");
		v_c05f[k]->get_linearization_coefficients(sub, range);
    //printf("qui\n");
	}
	else
		throw(std::logic_error("It's impossible to compute the subgradient because only diagonal linearization is allowed"));
  }
  sgb=nullptr;
	return(stp-strt);
  //printf("stop-start=%d\n", stp-strt);
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
  
   *f_log <<std::endl << "Stepsize =" << Stepsize_name << ", Deflection= " << Deflection_name<< std::endl;
   *f_log << std::endl <<"Scalls-Pariter="<< "{" << SCalls << "-" << ParIter << "} ";

   *f_log << " FiBar =";
    if( FiBar == Inf<VarValue>() )
     *f_log << "INF";
    else
      *f_log << -FiBar;

    *f_log << " ~ FiBest =";
    if( FiBest == Inf<VarValue>() )
     *f_log << "INF";
    else {
     *f_log << -FiBest;

     if( Lowerbound > - Inf<VarValue>() )
      *f_log << " ~ Lowerbound = "<< Lowerbound << " ~ Gap = " << ( FiBest - Lowerbound ) /
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
 
 *f_log << " FiLambda = ";  // print function value
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

 VLOG( 1 , " ~ Sigma1 = " << Sigma );

 }  // end( SubGradSolver::UpdateSigma() )

/*--------------------------------------------------------------------------*/

void SubGradSolver::UpdtLowerBound(void){
	//if the objective function is concave i take as lower bound its -upper bound 
	VarValue Lb= IsConvex ? f_Block->get_valid_lower_bound(false) : -(f_Block->get_valid_upper_bound(false));
		if(Lb> -Inf<VarValue>()){
      //printf("true lb");
			Lowerbound=Lb;
			TrueLB=true;
		}
	else {
		Lb= IsConvex ? f_Block->get_valid_lower_bound(true) : -(f_Block->get_valid_upper_bound(true));
    //printf("falselb\n");
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
 if(UBc!=0)
	tUbv.resize(num);
 else
	tUbv.resize(Ubv.size()+ strt);
 if(LBc!=0)
	tLbv.resize(num);
 else
	tLbv.resize(Lbv.size()+strt);

 tLbv.assign(num, -Inf<VarValue>());
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
  /*for(Index i=0; i<num; i++){
    //printf("Lbv[%d]=%f e Ubv[%d]=%f\n", i, Lbv[i], i , Ubv[i]);
  }*/
 //printf("Lbv.size()=%d", Lbv.size());
 /*std::vector <int>::size_type low=Lbv.size();
  std::cout <<" low=" << low << std::endl;
 std::vector <int>::size_type up=Ubv.size();
 std::cout << "up=" << up << std::endl ;*/
 }  // end( SubGradSolver::AllocateUBLB() )

/*--------------------------------------------------------------------------*/

void SubGradSolver::EvalGi( Index wFi ){
 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( wFi > NrFi )
  wFi = NrFi;//It must be one of the possible index


 Index* SGBse;
 Index SGBDim = GetGi( Gi.data() , SGBse , wFi?Inf<Index>():NrFi, NrFi );
  for(Index i=0; i<NumVar; i++){
    printf("In Evalgi Gi[%d]=%f\n", i , Gi[i]);
    //Gi[i]=gi[i];
  }
 // compute the linearization error of Gi[wFi] at Lambda1 - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //for(Index i=0; i<NumVar; i++)
  //printf("G[%d]=%f\n", i, Gi[i]);

 //VarValue upest=v_c05f[0]->get_upper_estimate();
 //printf("upest=%f\n", upest);
 //VarValue lowest=v_c05f[0]->get_lower_estimate();
 //printf("lowest=%f\n", lowest);
 Sigma=0;
 for(Index k=0; k<NrFi;k++)
 Sigma = Sigma + v_c05f[k]->get_linearization_constant();
 printf("Sigma in EvalGi=%f\n", Sigma);
 
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
printf("\n In ProjFsb:\n");
  std::vector <int>::size_type low=Lbv.size();
  std::vector <int>::size_type up=Ubv.size();
//   std::cout <<" low=" << low << std::endl;
//    std::cout <<" up=" << up << std::endl;
 if( ( low==0 ) && ( up==0 ) ){
   printf("caso low=up=0\n");
  return( LptHasChgd );
 }

 // project Lpt on the box constraints  - - - - - - - - - - - - - - - - - -

 if( (low>0) && (up>0) ) {
   printf("low>0 and up>0\n");
  // meanwhile, check that they are valid
  Index i = 0;
 for(; i<NumVar; i++)
  printf("Lbv[i]=%f\n", Lbv[i]);
  for(i=0; i<NumVar; i++)
  printf("Ubv[i]=%f\n", Ubv[i]);
  for( i=0; i < NumVar ; i++ ) {
   if( Ubv[ i ] - Lbv[ i ] < EpsLin ) {	
     printf("Caso Emptyset=true \n");
    EmptySet = true;
    break;
    }
   if( Lpt[ i ] > Ubv[ i ] ) {
     printf("Caso Lpt[i]>Ubv[i]\n");
    Lpt[ i ] = Ubv[ i ];
    LptHasChgd = true;
    break;
    }
   else
    if( Lpt[ i ] < Lbv[ i ] ) {
      printf("Caso Lpt[i]>Lbv[i]\n");
     Lpt[ i ] = Lbv[ i ];
     LptHasChgd = true;
     break;
     }
   }

  if( ! EmptySet ){
    printf("Caso emptyset\n");
   for(i=0 ; i < NumVar ; i++ ) {
    if( Ubv[ i ] - Lbv[ i ] < EpsLin ) {	
      printf("Caso Ubv[i]lLbv[i]<EpsLin\n");
     EmptySet = true;
     break;
     }

    if( Lpt[ i ] > Ubv[ i ] ){
     Lpt[ i ] = Ubv[ i ];
     printf("Lpt[i]>Ubv[i]\n");
    }
    else
     if( Lpt[ i ] < Lbv[ i ] ){
      Lpt[ i ] = Lbv[ i ];
      printf("Lpt[i]<Lbv[i]\n");
     }
    }
    }
  }
 else
  if( low>0 ) {
   Index i = 0;
   for( ; i < NumVar ; i++ )
    if( Lpt[ i ] < Lbv[ i ] ) {
      printf("low>0 Lpt[i]<Lbv[i]\n");
     Lpt[ i ] = Lbv[ i ];
     LptHasChgd = true;
     }
   /*for( i=0; i < NumVar ; i++ )
    if( Lpt[ i ] < Lbv[ i ] ){
      printf("low>0 Lpt[i]<Lbv[i]\n")
     Lpt[ i ] = Lbv[ i ];
    }*/
   }
  else {
   Index i = 0;
   for( i=0 ; i < NumVar ; i++ )
    if( Lpt[ i ] > Ubv[ i ] ) {
      printf("up>0 Lpt[i]>Ubv[i]\n");
     Lpt[ i ] = Ubv[ i ];
     LptHasChgd = true;
     }
   /*for( i=0; i < NumVar ; i++ )
    if( Lpt[ i ] > Ubv[ i ] )
     Lpt[ i ] = Ubv[ i ];
     */
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
	
 if( ( low==0 ) && ( up==0) )
  return;

 // project Gpt on the box constraints  - - - - - - - - - - - - - - - - - -

 if( (low>0) && (up>0) ) {
  for( Index j = strt ; j < stp ; j++ )
   if( ( ( Lpt[ j ] >= Ubv[ j ] ) && ( Gpt[ j ] <= 0 ) ) ||
       ( ( Lpt[ j ] <= Lbv[ j ] ) && ( Gpt[ j ] >= 0 ) ) )
    Gpt[ j ] = 0;
  }
 else{
  if( low>0 ) {
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
	for(Index i=0; i<n; i++){
		loc=loc+(A[i]*A[i]);
      //printf("loc=%f\n", loc);
  }

	return loc;
}

/*--------------------------------------------------------------------------*/

void SubGradSolver::InitStepsize( void ){
 step = 0;
 if (Stepsize_name=="FumeroTV")
 	stepsize=Stepsize::new_Stepsize("FumeroTV", this);
 if (Stepsize_name=="ColorTV")
 	stepsize=Stepsize::new_Stepsize("ColorTV", this);
  if (Stepsize_name=="Polyak")
  	stepsize=Stepsize::new_Stepsize("Polyak", this);
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

 for(i=0;i<NumVar;i++)
 //printf("dir[%d]=%f\n", i, dir[i]); 

 // if the pointer to Deflection object is not 0, invoke Format() in
 // Deflection, otherwise some conditions have to be taken into account:
 // i) the stepsize-restricted scheme must be used and ii) the projection
 // of the direction is useless
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if(Deflection_name=="Volume"){
  	deflection=Deflection::new_Deflection("Volume", this);
  	deflection->Format();
  }
   if(Deflection_name=="null"){
  	SGPar3 |= StepRST;//or bitwise
  	SGPar1 &= ~( PrjDIRM1 | PrjDIR );
  }
   if(Deflection_name=="PrimalDual"){
    printf("PrimalDual\n");
    deflection=dynamic_cast< PrimalDualClass *>(Stepsize::new_Stepsize("PrimalDual", this));
    printf("Dopo\n");
  	/*if(stepsize){
  	deflection=dynamic_cast< PrimalDual *>( stepsize );
  	if(!deflection){*/
  	deflection->Format();
    }
  /*	}
  	else
  		throw(std::invalid_argument("Deflection method not allowed"));
  }*/
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

VLOG( 1 , " ~ HFi = " << -FiHat );

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
 return( Solver->NrmDir);
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
 	return(Solver->v_c05f[0]->get_Lipschitz_constant());
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
 	return(Solver->v_c05f[0]->get_Lipschitz_constant());
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

PrimalDualClass::PrimalDualClassFactoryMap & PrimalDualClass::f_factory( void ){
	static PrimalDualClassFactoryMap s_factory;
	return( s_factory );
	}


/*---------------------------------------------------------------------------*/
/*--------------------------End File SubGradSolver.cpp-----------------------*/
/*---------------------------------------------------------------------------*/
