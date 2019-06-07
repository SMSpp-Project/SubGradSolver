/*
 * SubGradSolver.cpp
 *
 *  Created on: Apr 17, 2019
 *      Author: simone
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/


#include <SubGradSolver.h>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SubGradSolver to the Solver factory
SMSpp_insert_in_factory_cpp_0( SubGradSolver );


/*--------------------------------------------------------------------------*/
// define and initialize here the vector of int parameters names
const std::vector< std::string > SubGradSolver::int_pars_str =
			{ "intSGPar1" , "intSGPar3" , "intSGPar4","intSGPar5" };

// define and initialize here the vector of double parameters names
const std::vector< std::string > SubGradSolver::dbl_pars_str =
			{ "dblSGPar2" ,"dbltStar", "dblEpsLin","dblEInit","dblEFnalt","dblEDcrs","dblEStps" };



// define and initialize here the map for int parameters names
const std::map< std::string , SubGradSolver::idx_type > SubGradSolver::int_pars_map =
                   { { "intSGPar1"  , SubGradSolver::intSGPar1  } ,
		     { "intSGPar3" , SubGradSolver::intSGPar3 } ,
		     { "intSGPar4" , SubGradSolver::intSGPar4 } ,
		     { "intSGPar5" , SubGradSolver::intSGPar5 }};

// define and initialize here the map for double parameters names
const std::map< std::string , SubGradSolver::idx_type > SubGradSolver::dbl_pars_map =
                   { { "dblSGPar2" , SubGradSolver::dblSGPar2 } ,
		     { "dbltStar" , SubGradSolver::dbltStar } ,
			 { "dblEpsLin" , SubGradSolver::dblEpsLin } ,
			 { "dblEInit" , SubGradSolver::dblEInit } ,
			 { "dblEFnalt" , SubGradSolver::dblEFnalt } ,
			 { "dblEDcrs" , SubGradSolver::dblEDcrs } ,
			 { "dblEStps" , SubGradSolver::dblEStps } };


// define and initialize here the default int parameters
const std::vector<int> SubGradSolver::dflt_int_par =
        {    10 ,  // intSGPar1
			100 ,  // intSGPar3
			  0 ,  // intSGPar4
			  0   // intSGPar5
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

/*--------------------------------------------------------------------------*/


#if ! SubGrad_HANDLES_CONSTRAINTS
 #define MaxName 3
 #define MaxNCnst 0
#endif


/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static const double Multplier = 1.0;

static const unsigned char PrjSG =  1;   // project of Gi[wFi] = \f$g_i\f$
static const unsigned char PrjDIRM1 = 2; // project of DirM1 = \f$d_{i-1}\f$
static const unsigned char PrjDIR = 4;   // project of DirM1 = \f$d_i\f$

static const unsigned char SafeRULE = 1;  // safe rule
static const unsigned char StepRST = 2;   // stepsize-restricted

static const unsigned char RstAlg =  1;  // don't reset algorithmic parameters
static const unsigned char RstCrr =  2;  // don't reset current point
static const unsigned char RstSbg =  4;  // don't reset subgradients
static const unsigned char RstCnt =  8;  // don't reset constraints
static const unsigned char RstFiV = 16;  // don't reset FiVals

static const bool first = true;          // indicate the position of the
static const bool second = false;        // direction

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

int myrandom( int i ) { return( rand() % i ); }

/*--------------------------------------------------------------------------*/




SubGradSolver::SubGradSolver() {
	 // initialize algorithmic parameters - - - - - - - - - - - - - - - - - - - -

	KSGPar1=get_int_par(int_par_str2idx("intSGPar1"));
	KSGPar2=get_dbl_par(dbl_par_str2idx("dblSGPar2"));
	KSGPar3=get_int_par(int_par_str2idx("intSGPar3"));
	KSGPar4=get_int_par(int_par_str2idx("intSGPar4"));
	KSGPar5=get_int_par(int_par_str2idx("intSGPar5"));
	srand( KSGPar5 );





	 //if( SGPar2 < 0 )

	//	 this->set_par(this->dbl_par_str2idx("SGPar2"),double(0));


	 // some initializations - - - - - - - - - - - - - - - - - - - - - - - - - - -

//stepsize = nullptr;
	 //deflection = nullptr;
	 //Lambda = LambdaHat = LambdaBest = nullptr;
	 //Gi = dir = nullptr;
	//lb =- Inf<double>();
	//ub = double (0);

//	 #if SubGrad_HANDLES_CONSTRAINTS
//	  MaxNCnst = 0;
//	  Q2KNP = nullptr;
//	  CnstBeg = nullptr;
//	  CnstVol = nullptr;
//	  CnstNxt = nullptr;
//	 #endif

	 KpBstL = false;
	 LHasChgd = true;  // LHasChgd will be again set to true by FormLambda1



}

SubGradSolver::~SubGradSolver() {
	// TODO Auto-generated destructor stub
}



void SubGradSolver::set_Block(Block* block) {


	 Solver::set_Block( block );  // attach to the new Block

	 /* Two types of block can be handled by the BundleSolver:

	     1. Only one single non-smooth function
	     2. a sum of non-smooth functions

	    The algorithm here developed aims at solving non-constrained non-smooth
	    optimization. The block can have box constraints at the most.
	    It is expected the block to have in the first case a FRealObjective
	    whose the function is a C05Function one and having no children, while in
	    the second case a FRealObjective one whose the function is a LinearFunction
	    and having as many sub-blocks as the number of components. In the latter case,
	    each sub-block must not contain any Variable or Constraint. */

	 if( f_Block->get_nested_Blocks().empty()  ) {

	   // the objective function of the block must be a C05Function  - - - - - - - -
	   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	   auto obj = boost::any_cast<FRealObjective *>( f_Block->get_objective() );
	   if( obj == nullptr )
	    throw( std::logic_error( "the objective is not a real function" ) );

	   auto c05f = dynamic_cast<C05Function *>( (obj)->get_function() );
	   if( c05f == nullptr )
	    throw( std::logic_error( "the objective is not a C05Function" ) );
	   v_c05f.push_back( c05f );

	   lf = nullptr;

	   }
	  else {

	   // the objective function of each block must be a LinearFunction - - - - - -
	   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	   auto obj = boost::any_cast<FRealObjective *>( f_Block->get_objective() );
	   if( obj == nullptr )
	    throw( std::logic_error( "the objective is not a real function" ) );

	   lf = dynamic_cast<LinearFunction *>( (obj)->get_function() );
	   if( lf == nullptr )
	    throw( std::logic_error( "the objective is not a LinearFunction" ) );

	   auto sb = f_Block->get_nested_Blocks();
	   v_c05f.resize( sb.size() );

	   for( Variable::Index i = 0 ; i < sb.size() ; ++i ) { // for each sub-block

	    // the objective function of each sub-block must be a C05Function - - - - -
	    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	    auto obj = boost::any_cast<FRealObjective *>( sb[ i ]->get_objective() );
	    if( obj == nullptr )
	     throw( std::logic_error( "the objective is not a real function" ) );

	    auto c05f = dynamic_cast<C05Function *>( (obj)->get_function() );
	    if( c05f == nullptr )
	 	throw( std::logic_error( "the objective is not a C05Function" ) );
	    v_c05f[ i ] = c05f;

	    // nephew are not allowed - - - - - - - - - - - - - - - - - - - - - - - - -
	    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	    if( sb[ i ]->get_nested_Blocks().size() )
	 	throw( std::logic_error( "nephew are not allowed" ) );

	    // Variable of Sub-Block are not expected, neither the Constraint - - - - -
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
	   } // end decomposed case - - - - - - - - - - - - - - - - - - - - - - - - - -

	  // if some Variable are present, they are of the ColVariable type - - - - - -
	  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  NumVar = 0;  // count the number of Variable
	  auto v_s_Variable = f_Block->get_static_variables();
	  for( auto & el : v_s_Variable ) {
	  if( un_any_thing_0( ColVariable , el , [ & ]{ ++NumVar; } ) )
	   break;
	  if( un_any_thing_1( ColVariable , el , [ & ]{ NumVar += var.size(); } ) )
	   break;
	  if( un_any_thing_K( ColVariable , el , [ & ]{ NumVar += var.size(); } ) )
	   break;
	  throw( std::logic_error( "some static Variable is not a ColVariable" ) );
	  }

	  // construct the vocabulary for Variable and sort it  - - - - - - - - - - - -
	  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  LamVcblr.resize( NumVar );
	  int count = 0;
	  for( auto & el : v_s_Variable )
	   un_any_static( el , [ & ]( ColVariable & static_var ){ LamVcblr[count++]= &static_var; } ,
	 		  un_any_type<ColVariable>() );

	  std::sort( LamVcblr.begin() , LamVcblr.end() );

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



	  //setfi

	  //MemDealloc();  // deallocate memory


	  MaxNumVar = NumVar//Oracle->GetMaxNumVar();


	//nella mail dice che a simone non interessano
	//if( NrFi > 1 )
	     //for( Variable::Index k = 0 ; k++ < NrFi ; )
	      //if( Oracle->GetBNC( k ) )
	       //throw( NDOException(
	  	//	   "SubGrad::SetFiOracle: easy component not allowed yet" ) );

	    // allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		std::vector<double> vd0(MaxNumVar,0);
	    LambdaBar.swap(vd0); //= new std::vector<double>(MaxNumVar) LMNum[ MaxNumVar ];                // stability center



	   // VectAssign( LambdaBar , LMNum( 0 ) , MaxNumVar );  // default = 0

	    std::vector<double> vd0(MaxNumVar,0);
	    Lambda.swap(vd0);
	    //Lambda = new LMNum[ MaxNumVar ];                   // tentative point
	    //VectAssign( Lambda , LMNum( 0 ) , MaxNumVar );     // default = 0

	    if( KpBstL ){
	    	std::vector<double> vd0(MaxNumVar,0);
	    	LambdaBest.swap(vd0);            // best point so far
	    }
	    else
	     LambdaBest.clear();

	    LowerBound = -Inf<double>();               // lower bound on the function
	    TrueLB = false;

	    FiBest = FiLambda = Inf<double>();         // Fi( Lambda ) is not known

	    // allocate memory for the subgrads and the direction - - - - - - - - - - -
	    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	    // The direction will be allocated after either Deflection will have been
	    // defined or changed [see InitDeflection()]. Instead, the subgradient is
	    // considered separately from the deflection and it will keep on having the
	    // same values after any deflection update.
	    std::vector<double> g0(MaxNumVar,0);
	    std::vector<double> d0(MaxNumVar,0);
	    Gi.swap(g0);
	    dir.swap(d0);

	    NrmGi = 0;

	    // allocate memory for box constraints- - - - - - - - - - - - - - - - - - -
	    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



	    ///////to DOOOOOO
//	    AllocateUBLB( 0 , NumVar );
//
//	    // allocate memory for disjoint simplex constraints - - - - - - - - - - - -
//	    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//	    #if ! SubGrad_HANDLES_CONSTRAINTS
//	     if( Oracle->MaxNConst() )
//	      throw( NDOException( "SubGrad::SetFiOracle: non-box constraints allowed"
//	  			  ) );
//	    #else
//	     MaxNCnst = Oracle->MaxNConst();
//
//	     if( MaxNCnst ) {
//	      CnstBeg = new Index[ MaxNCnst ];
//	      CnstNxt = new SIndex[ MaxNumVar ];
//	      CnstVol = new HpNum[ MaxNCnst ];
//	      VectAssign( CnstBeg , Index( Inf<Index>() ) , MaxNCnst );
//	      VectAssign( CnstNxt , SIndex(-1) , MaxNCnst );
//	      MaxName = MaxNCnst + 3;
//	      }
//	     else
//	      MaxName = 3;
//	    #endif
//
//	    ReSetAlg( RstCrr | RstSbg | RstCnt );
//
//	    // tell the oracle about the NDOSolver and its settings - - - - - - - - - -
//	    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//	    Oracle->SetNDOSolver( this );
//	    Oracle->SetMaxName( MaxName );
//	    Oracle->SetPrecision( EInit );
//
//	    // warning: the following things can only be done *after* that
//	    // Oracle->SetMaxName() has been invoked, because they use methods of the
//	    // oracle which depends on knowledge of the MaxName to work properly
//	    // insert the constant subgradient of the 0-th component - - - - - - - - - -
//
//	    MultBse = new Index[ 2 ];
//	    MultBse[ 1 ] = Inf<Index>();
//
//	    // asking the projection of both DirM1 and and Gi implies the automatic
//	    // projection of dir- - - - - - - - - - - - - - - - - - - - - - - - - - - -
//	    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//	    if( ( SGPar1 & PrjDIRM1 ) && ( SGPar1 & PrjSG ) )
//	     SGPar1 &= ~PrjDIR;
//
//	    // checking for 0-th component - - - - - - - - - - - - - - - - - - - - - - -
//	    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//	    cIndex_Set SG0Bse;
//	    cIndex SG0BDm = Oracle->GetGi( Gi , SG0Bse , MaxName );
//
//	    ZeroComp = SG0BDm > 0;
//	    SubGrad::SetPar( SubGrad::kSGPar2 , SGPar2 );
//	    // initialization of data structures for incremental done there inside
//
//	    // initialize HatLambda  - - - - - - - - - - - - - - - - - - - - - - - - - -
//	    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//	    if( SGPar4 ) {
//	     LambdaHat = new LMNum[ MaxNumVar ];
//	     VectAssign( LambdaHat , LMNum( 0 ) , MaxNumVar );
//	     }
//	    else
//	     LambdaHat = nullptr;
//
//	    FiHat = 0;     // Fi( HatLambda ) is not known
//
//	    // initialize both stepsize and deflection solvers - - - - - - - - - - - - -
//	    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//	    if( stepsize )
//	     InitStepsize();
//
//	    InitDeflection();
//
//	    // other initializations- - - - - - - - - - - - - - - - - - - - - - - - - -
//
//	    fs = FiOracle::kFiNorm;
//	    EmptySet = false;       // the feasible set is not empty
//
//	    } // end if( Fi )
//	   else {  // setting a null FiOracle
//	    if( deflection )
//	     deflection->Format();
//
//	    if( stepsize )
//	     stepsize->Format();
//	    }




}

void SubGradSolver::set_log(std::ostream* log_stream) {
	 Solver::set_log(log_stream);

	 *f_log << std::endl << "SubGrad Solver: Vars = " << "NumVar" << " ("
	 	   << MaxNumVar << ") ~  NItIncr = " << NItIncr <<  " ~ EpsLin = "
	 	   << EpsLin << " ~ FiEps = " << "EInit" << " ~ t* = " << tStar;

	    *f_log << std::endl << " SGParameters = ( " << KSGPar1 << " ~ " << KSGPar2
	 	   << " ~ " << KSGPar3 << " ~ " << ( KSGPar4? "true" : "false")
	 	   <<  " ~ "  << KSGPar5 << " )" << std::endl;


 }




int SubGradSolver::compute(bool changedvars) {

	// initializations- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 Result = kOK;
	 SCalls++;
//	 if( NDOt )
//	  NDOt->Start(); //timer

	 Variable::Index comp; // a component

	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 do {

	  // function log: print the information relative to the stability center
	  // and the best point - - - - - - - - - - - - - - - - - - - - - - - - - - -

		 Log1();


	 // check for running time - - - - - - - - - - - - - - - - - - - - - - - - -
	   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



		 // To compute the next point Lambda1, select the direction Dir and the
		   // stepsize Nu using, respectively, Deflection and Stepsize objects
		   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

		   if( FiBar < Inf<double>() ) {

			   // if the function value is unknown  skip the direction and stepsize
			      // computation and go straight to FiAndGi()

			      // direction log: print the information relative to DIRM1  - - - - - - - -
			      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			      Log2();

			      // FormD call: a new direction is provided, but Dir is not updated
			      // yet because AddVariables() and RemoveVariables() need to work
			      // with the prior vector Dir - - - - - - - - - - - - - - - - - - - - - - -
			      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			      FormD( );

			         // check the status of the FiOracle and take the action, however
			         // because we do not succeed in handling the variables generation
			         // combined with incremental approach, GetFiStatus is invoked only during
			         // outer iterations  - - - - - - - - - - - - - - - - - - - - - - - - - - -
			         // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


			      if( !InnIter )
			          //fs = Oracle->GetFiStatus();

			         if( fs == kFiStop ) {
			          VLOG( 1 , " ~ FiOracle:STOP" << std::endl );
			          Result = kStopIter;
			          break;
			          }

			         if( fs == kFiError ) {
			          VLOG( 1 , " ~ Error in the FiOracle" << std::endl );
			          Result = kError;
			          break;
			          }

			         // if something is changed go back to FormD(): the deflection coefficient
			         // Alfa, and in turn Dir, must be again computed, and, in the
			         // deflection-restricted approach, even the stepsize Nu- - - - - - - - - -
			         // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			         if( fs == kFiChgd ) {
			          VLOG( 1 , " ~ Fi changed: loop" << std::endl );
			          continue;
			          }

			         // save the search direction Dir: note that the previous DIR is not saved
			         // as DirM1, but it is merely discharged; only and one vector is kept
			         // for storing the direction- - - - - - - - - - - - - - - - - - - - - - - -
			         //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			         SaveDir();


			         // if the Stepsize object needs the scalar product between DirM1 and Gi,
			            // then the direction vector is kept- - - - - - - - - - - - - - - - - - - -

			            if( stepsize->NeedsdkM1Gk() && ( ! deflection ) )
			             dir.swap(Gi);


			            // direction log: print the information relative to DIR   - - - - - - - - -
			            //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			            Log2();

			            // a stop condition regarding consecutive "short" stepsizes - - - - - - - -
			            //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			            if( step <= 1e-8 * tStar ) {
			             CSmallStep++;
			             if( CSmallStep > ( ( NItIncr > 0 )? 100 * NrFi : 100 ) ) {
			              Result = kStopIter;
			              break;
			              }
			             }
			            else
			             CSmallStep = 0;

			            }  // end( search direction computation ) - - - - - - - - - - - - - - - - -
			               // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           // check for optimality  - - - - - - - - - - - - - - - - - - - - - - - - - -
			           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           if( ( ! InnIter ) /*&& IsOptimal()*/ )
			            break;

			           // calculate the tentative point Lambda1 (not projected yet)   - - - - - - -
			           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           FormLambda1();

			           VLOG( 1 , std::endl << "           " << " |Lambda1| = "
			         	         << std::inner_product(Lambda.begin(),Lambda.end(),Lambda.begin(),0));

			           // a real iteration (iterations where Fi() is not evaluated count as well) -
			           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           ParIter++;

			           // is an inner iteration?- - - - - - - - - - - - - - - - - - - - - - - - - -
			           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           InnIter = NItIncr? GiEvaltns % ( NItIncr + 1 ): 0;
			           LHasProj = false;

			           // pick up a component - - - - - - - - - - - - - - - - - - - - - - - - - - -
			           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           Variable::Index comp;
			           if( InnIter ) {  // if it is the first component, re-shuffle the vector
			            comp = ( ( GiEvaltns % ( NItIncr + 1 ) ) - 1 ) %
			                     ( ZeroComp? NrFi + 1 : NrFi );
			            if( ! comp )
			             std::random_shuffle( Seq.begin() , Seq.end() , myrandom );
			            comp = Seq.at( comp );
			            }
			           else
			            comp = Inf<Variable::Index>();

			           // project Lambda1 onto the feasible set - - - - - - - - - - - - - - - - - -
			           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           for( ;; ) {
			            if( ProjFsb( Lambda , EmptySet ) ) {  // go to the feasible point - - - - -
			             LHasProj = LHasChgd = true;
			             VLOG( 1 , " ~ |Lambda1| = " << std::inner_product(Lambda.begin(),Lambda.end(),Lambda.begin(),0));
			             if( EmptySet )
			              break;
			             }

			            // never evaluate the objective function at the same point - - - - - - - -
			            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			            if( ! LHasChgd )
			             throw( std::logic_error( "SubGrad::Solve: this should not happen" ) );

			            // compute Fi( Lambda1 ) and a subgradient / constraint Gi[wFi]- - - - - -
			            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			            FiAndGi( comp );

			            // if the Lambda1 is infeasible, go back to the projection - - - - - - - -
			            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			            if( FiLambda < Inf<double>() )
			             break;

			            } // end function evaluation - - - - - - - - - - - - - - - - - - - - - - -
			              // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           // project the subgradient Gi[wFi] over the tangent cone at tentative point
			           // Lambda1- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			           //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           if( KSGPar1 & PrjSG )
			            ProjTng( Gi , Lambda );

			           NrmGi = std::inner_product(Gi.begin(),Gi.end(),Gi.begin(),0);

			           // compute the scalar product between Dir and Gi[wFi] at Lambda1- - - - - -
			           //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           if( stepsize->NeedsdkM1Gk() ) {
			            dM1Gk = std::inner_product(dir.begin(),dir.end(),Gi.begin(),0);
			            dM1GkDone = true;
			            }
			           else
			            dM1GkDone = false;

			           // print the logs regarding Fi and Gi[wFi] at Lambda1 - - - - - - - - - - -
			           //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           Log3();
			           Log4();

			           // check for feasibility- - - - - - - - - - - - - - - - - - - - - - - - - -
			           //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           if( EmptySet ) {
			            Result = kInfeasible;// kUnfsbl;
			            break;
			            }

			           if( FiLambda == -Inf<double>() ) {
			            Result = kUnbounded;//kUnbndd;
			            break;
			            }

			           // check whether or not the Lower Bound has changed - - - - - - - - - - - -
			           //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           //UpdtLowerBound();
			           VLOG2( 1 , LowerBound > - Inf<double>() , " ~ UB = " << -LowerBound );

			           // if FiBar (the function at LambdaBar) is unknown, or no deflection is
			           // used, then move to Lambda1 - - - - - - - - - - - - - - - - - - - - - - -
			           //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           if( ( FiBar == Inf<double>() ) || ( ! deflection ) ) {
			            GotoLambda1();

			            if( deflection ) {
			             double Mlt = 0;
			             Variable::Index NmSt = MaxNCnst;


			             C05Function::LinearCombination lc({std::make_pair(0,Mlt)});




					   for(auto c05 : v_c05f)
						   c05->store_combination_of_linearizations(lc , DirPos == first ? 2 : 1);

					   if(lf)
						   lf->store_combination_of_linearizations(lc , DirPos == first ? 2 : 1);




			            // Oracle->Aggregate( &Mlt , &NmSt , 1 , MaxNCnst + 2 );
			             }

			            DoSS = true;

			            VLOG( 1 , std::endl );
			            continue;
			            }

			           // Null Step (NS) vs Serious Step (SS)- - - - - - - - - - - - - - - - - - -
			           //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           if( ( DoSS = deflection->DoSS() ) ) {  // SS- - - - - - - - - - - - - - - -

			            VLOG2( 1 , deflection->Delta() < Inf<double>() , std::endl << " SS["
			         	      << CSSCntr );

			            GotoLambda1();

			            CSSCntr++;
			            CNSCntr = 0;
			            }
			           else {  // NS- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			            VLOG2( 1 , deflection->Delta() > - Inf<double>() , std::endl << " NS["
			         	      << CNSCntr );

			            UpdateSigma();

			            CNSCntr++;
			            CSSCntr = 0;

			            }  // end SS/NS - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			               // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           VLOG( 1 , "]: DFi (" << FiBar - FiLambda << ") < m1 * Dv ("
			         	     << deflection->Delta() << ")" );
			           VLOG( 1 , std::endl );

			           // check for Lower Bound stopping condition - - - - - - - - - - - - - - - -
			           //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			           if( TrueLB )
			            if( FiBest - EpsLin * std::abs(FiBest ) <= LowerBound )
			             break;

			           if( FiBest <= LowerBound  * ( 1.0 - ( LowerBound > 0.0 ? EpsLin : - EpsLin )
			         				) ) {
			            Result = kUnbounded;//kUnbndd;
			            break;
			            }
			           } while( ( ! MaxIter ) || ( ParIter < MaxIter ) );

			          // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			          // main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
			          // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

			          if( MaxIter && ( ParIter >= MaxIter ) && ( ! Result ) )
			           Result = kStopIter;//kStpIter;

//			          if( NDOt )
//			           NDOt->Stop();

			          return( Result );


}





//}

void SubGradSolver::get_var_solution(Configuration* solc) {
}

void SubGradSolver::get_dual_solution(Configuration* solc) {

	for( Variable::Index i = 0 ; i < zA.size() ; ++i ) {
	  if( zA[i].second.empty() )
	   throw( std::invalid_argument( "the combination is not present" ) );
	 // v_c05f[ i ]->set_important_linearization( zA[i].first );
	  }

}




//
//
//void SubGradSolver::SetDeflection(Deflection* Vol) {
//}
//



int SubGradSolver::get_int_par(const idx_type par) const {


	switch( par ) {
	  case( intMaxIter ):
	   return( MaxIter );
	   break;
	  case( intMaxSol ):
	   return( MaxSol );
	   break;
	  case( intLogVerb ):
	   return( LogVerb );
	   break;
	  case( intSGPar1 ):
	   return( KSGPar1 );
	   break;
	  case( intSGPar3 ):
	   return( KSGPar3 );
	   break;
	  case( intSGPar4 ):
	   return( KSGPar4 );
	   break;
	  case( intSGPar5 ):
	   return( KSGPar5 );
	   break;
	  default:
	     return( get_dflt_int_par( par ) );
	}

}

double SubGradSolver::get_dbl_par(const idx_type par) const {

	switch( par ) {
	  case( dblMaxTime ):
	   return( MaxTime );
	  case( dblRelAcc ):
	   return( RelAcc );
	  case( dblAbsAcc ):
	   return( AbsAcc );
	  case( dblUpCutOff ):
	   return( UpCutOff );
	  case( dblLwCutOff ):
	   return( LwCutOff );
	  case( dblRAccSol ):
	   return( RAccSol );
	  case( dblAAccSol ):
	   return( AAccSol );
	  case( dblFAccSol  ):
	   return( FAccSol );
	  case( dbltStar ):
	   return( tStar );
	  case( dblEInit ):
	   return( EInit );
	  case( dblEFnalt ):
	   return( EFnal );
	  case( dblEDcrs ):
	   return( EDcrs );
	  case(dblSGPar2):
		return(KSGPar2);
	  default:
	    return( get_dflt_dbl_par( par ) );
	}


}


void SubGradSolver::set_par(const idx_type par, const int value) {

	switch( par ) {
	  case( intMaxIter ):
	    MaxIter =value;
	   break;
	  case( intMaxSol ):
	    MaxSol =value;
	   break;
	  case( intLogVerb ):
	   LogVerb =value;
	   break;
	  case( intSGPar1 ):
	   KSGPar1 =value;
	   break;
	  case( intSGPar3 ):
	   KSGPar3 =value;
	   break;
	  case( intSGPar4 ):
	   KSGPar4 =value;
	   break;
	  case( intSGPar5 ):
	   KSGPar5 =value;
	   break;
	}
}

void SubGradSolver::set_par(const idx_type par, const double value) {
	switch( par ) {
	  case( dblMaxTime ):
	    MaxTime =value;
		break;
	  case( dblRelAcc ):
	    RelAcc =value;
		break;
	  case( dblAbsAcc ):
	    AbsAcc =value;
	  	break;
	  case( dblUpCutOff ):
	    UpCutOff =value;
	    break;
	  case( dblLwCutOff ):
	   LwCutOff =value;
	    break;
	  case( dblRAccSol ):
	    RAccSol =value;
	    break;
	  case( dblAAccSol ):
	    AAccSol =value;
	    break;
	  case( dblFAccSol  ):
	    FAccSol =value;
	    break;
	  case( dbltStar ):
	    tStar =value;
	    break;
	  case( dblEInit ):
	    EInit = value;
	    break;
	  case( dblEFnalt ):
	    EFnal =value;
	    break;
	  case( dblEDcrs ):
	    EDcrs =value;
	    break;
	  case(dblSGPar2):
		KSGPar2 =value;
	    break;
	}
}














ComputeConfig* SubGradSolver::get_ComputeConfig(bool all,
		ComputeConfig* ocfg) const {
	auto tc=ThinComputeInterface::get_ComputeConfig(all,ocfg);

	auto sc = tc->new_Configuration("SimpleConfig_p_p");

	auto s= dynamic_cast<SimpleConfiguration<std::pair< Configuration * , Configuration * > >*>(sc);


	s->f_value.first=this->getConfstepsize();
	s->f_value.second=this->getConfdeflection();


	tc->f_extra_Configuration=s;




	return tc;


}




void SubGradSolver::set_ComputeConfig(ComputeConfig* scfg) {

	ThinComputeInterface::set_ComputeConfig(scfg);



}






  Configuration* SubGradSolver::getConfdeflection() const {
	return confdeflection;
}

void SubGradSolver::setConfdeflection(Configuration*& confdeflection)  {
	this->confdeflection = confdeflection;
}

  Configuration* SubGradSolver::getConfstepsize()  const{
	return confstepsize;
}



void SubGradSolver::SetStepsize(Stepsize* STP) {

//	if( stepsize ) {
//	  // if the Stepsize and Deflection objects are actually the same, ensure
//	  // it is not deleted twice
//	  auto StpDef = dynamic_cast<Deflection *>( stepsize );
//	  delete stepsize;
//	  if( ! StpDef )
//	   deflection = nullptr;
//	  }
//
//	 stepsize = STP;
//
//	 if( stepsize && Oracle )
//	  InitStepsize();

}

void SubGradSolver::SetDeflection(Deflection* Vol) {
	if( deflection ) {
	  // if the Stepsize and Deflection objects are actually the same, ensure
	  // it is not deleted twice
	  auto StpDef = dynamic_cast<Stepsize *>( deflection );
	  delete deflection;
	  if( ! StpDef )
	   stepsize = nullptr;
	  }

	 deflection = Vol;

	 if( ! v_c05f.empty() )
	  InitDeflection();
	 else
	  alpha = 1;
}

void SubGradSolver::SetFiOracle(C05Function* Fi) {
}

void SubGradSolver::SetLambda(std::vector<double> tLambda) {
}

void SubGradSolver::KeepBestLambda(const bool KBL) {
}

void SubGradSolver::setConfstepsize( Configuration*& confstepsize) {
	this->confstepsize = confstepsize;
}




Function::Vec_FunctionValue SubGradSolver::ReadBestSol(const std::vector<Variable::Index> &I , Variable::Index &D ) {

	//I = nullptr;
	D = NumVar;
	return( KpBstL ? LambdaBest : LambdaBar );

}

Function::FunctionValue SubGradSolver::ReadBestFiVal(Variable::Index wFi) {

	if( wFi!= Inf<double>() )
	  std::logic_error( "SubGrad::ReadBestFiVal: wFi != INF" );

	 return( FiBest );

}

Function::Vec_FunctionValue SubGradSolver::ReadSol(
		const std::vector<Variable::Index> &I , Variable::Index &D) {
	//I = nullptr;
	 D = NumVar;
	 return( LambdaBar );

}

Function::FunctionValue SubGradSolver::ReadFiVal( Variable::Index wFi ) {

	if( wFi!= Inf<Variable::Index>() )
	  std::logic_error( "SubGrad::ReadFiVal: wFi != INF"  );

	 return( FiBar );


}

Function::FunctionValue SubGradSolver::ReadHatFiVal(void) {

	if( intSGPar4 )
	  return( FiHat );
	 else
	  return( Inf<double>() );

}



const Function::FunctionValue* SubGradSolver::ReadMult( Variable::Index &I , Variable::Index &D , Variable::Index wFi ) {

//	if( wFi!= Inf<Variable::Index>() )
//	  std::logic_error( "SubGrad::ReadBestFiVal: wFi != INF" );
//
//	 if( FiBar < Inf<double>() ) {
//	  if( deflection )
//	    MultBse= ( DirPos == first ) ? ( MaxNCnst + 2 ) : ( MaxNCnst + 1 );
//	  else
//	   MultBse = MaxNCnst;
//
//	  I = MultBse;
//	  D = 1;
//	  return( & Multplier );
//	  }
//	 else {
//	 // I = nullptr;
//	  D = 0;
//	  return( 0 );
//	  }

}

const  Function::FunctionValue  SubGradSolver::ReadLBMult(const Variable::Index wFi) {

	if( wFi!= Inf<Variable::Index>() )
	   std::logic_error( "SubGrad::ReadBestFiVal: wFi != INF" );

	 return( 0 );

}

Variable::Index SubGradSolver::NrIter(void) const {
	return( ParIter );

}

const double SubGradSolver::ReadFVal( void){

	return FiLambda;

}

void SubGradSolver::AllocateUBLB(const Variable::Index strt,
		const Variable::Index num) {

//to do
//	std::vector<double> newp(num);
//	std::vector<double> tub = ub ? ub + strt : newp;
//	 LMRow tlb = lb.empty() ? lb + strt : new LMNum[ num ];
//	 VectAssign( tlb , - Inf<LMNum>() , num );
//
//	 bool LBc = lb;
//	 bool UBc = ub;
//
//	 Index i = 0;
//	 if( ( ! LBc ) && ( ! UBc ) )
//	  for( ; i < num ; ++i ) {
//	   if( ! Oracle->GetUC( strt + i ) ) {
//	    tlb[ i ] = 0;
//	    LBc = true;
//	    }
//
//	   if( ( tub[ i ] = Oracle->GetUB( strt + i ) ) < Inf<LMNum>() )
//	    UBc = true;
//
//	   if( LBc || UBc )
//	    break;
//	   }
//
//	 if( LBc && ( ! UBc ) )
//	  for( ; i < num ; ++i ) {
//	   if( ! Oracle->GetUC( strt + i ) )
//	    tlb[ i ] = 0;
//
//	   if( ( tub[ i ] = Oracle->GetUB( strt + i ) ) < Inf<LMNum>() )
//	    UBc = true;
//
//	   if( UBc )
//	    break;
//	   }
//
//	 if( ( ! LBc ) && UBc )
//	  for( ; i < num ; ++i ) {
//	   if( ! Oracle->GetUC( strt + i ) ) {
//	    tlb[ i ] = 0;
//	    LBc = true;
//	    }
//
//	   tub[ i ] = Oracle->GetUB( strt + i );
//
//	   if( LBc )
//	    break;
//	   }
//
//	 for( ; i < NumVar ; i++ ) {
//	  if( ! Oracle->GetUC( strt + i ) )
//	   tlb[ i ] = 0;
//	  tub[ i ] = Oracle->GetUB( strt + i );
//	  }
//
//	 if( ( ! lb ) && LBc ) {
//	  if( ( strt == 0 ) && ( num == MaxNumVar ) )
//	   lb = tlb;
//	  else {
//	   lb = new LMNum[ MaxNumVar ];
//	   cLMNum dflt = - Inf<LMNum>();
//	   VectAssign( lb , dflt , strt );
//	   VectAssign( lb + strt , tlb , num );
//	   delete[] tlb;
//	   }
//	  }
//
//	 if( ( ! ub ) && UBc ) {
//	  if( ( strt == 0 ) && ( num == MaxNumVar ) )
//	   ub = tub;
//	  else {
//	   ub = new LMNum[ MaxNumVar ];
//	   cLMNum dflt = Inf<LMNum>();
//	   VectAssign( ub , dflt , strt );
//	   VectAssign( ub + strt , tub , num );
//	   delete[] tub;
//	   }
//	  }


}

double SubGradSolver::lbnd(const Variable::Index i) {
	//return( lb ? lb[ i ] : - Inf<double>() );
}

double SubGradSolver::ubnd(const Variable::Index i) {
	//return( ub ? ub[ i ] : Inf<double>() );
}

void SubGradSolver::EvalGi(Variable::Index wFi) {

	 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( wFi > NrFi )
	  wFi = Inf<Variable::Index>();

//?????????????????????
//	 if( wFi && ! Oracle->get_linearization_constant( wFi ))// NewGi( wFi ) )
//	  throw( std::logic_error(
//			    "SubGradSolver::EvalGi: the subgradient cannot be computed" ) );

	 // fetch the item (subgradient) from the Oracle - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 //const std::vector<Variable::Index> SGBse;
	 //Variable::Index SGBDim = Oracle->get_linearization_coefficients() GetGi( Gi , SGBse , wFi? Inf<Variable::Index>() : MaxName );
	 std::vector<double> pg(NumVar);
	 	  for(auto c05 : v_c05f){
	 		  c05->get_linearization_coefficients(pg.data());
	 		  std::transform (Gi.begin(), Gi.end(), pg.begin(), Gi.begin(), std::plus<double>());


	 	  }
	 	  // if the case of when the 0-th component is a linear function
	 	  if(lf){
	 		  lf->get_linearization_coefficients(pg.data());
	 		  std::transform (Gi.begin(), Gi.end(), pg.begin(), Gi.begin(), std::plus<double>());

	 	  }



	 // densify the item  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 //if( SGBse )
	  //Densify( Gi , SGBse , SGBDim , NumVar );

	 // compute the linearization error of Gi[wFi] at Lambda1 - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 //Sigma = Oracle->get_value() ;//GetVal();

	 for(int i=0; i<v_c05f.size();i++){
		  v_c05f[i]->compute(true);
		  Sigma+=v_c05f[i]->get_value();
	  }

	 if( KSGPar4 )
	  SigmaHat = Sigma;

	 // write the item  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( wFi && ( ! InnIter ) )
		 for(auto a : v_c05f)
			 a->store_linearization(wFi);// SetGiName( MaxNCnst );

	 GiEvaltns++;

}

bool SubGradSolver::ProjFsb(std::vector<double> Lpt, bool& EmptySet) {

	bool LptHasChgd = false;  // it is true, if Lpt changes
	 EmptySet = false;         // it is true, if the feasible set is empty

	 // the case with knapsack constraints is treated apart - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 #if SubGrad_HANDLES_CONSTRAINTS

	 #endif

	 if( ( ! lb.empty() ) && ( ! ub.empty() ) )
	  return( LptHasChgd );

	 // project Lpt on the box constraints  - - - - - - - - - - - - - - - - - -

	 if( lb.empty() && ub.empty() ) {
	  // meanwhile, check that they are valid

	  Variable::Index i = 0;
	  for( ; i < NumVar ; i++ ) {
	   if( ub[ i ] - lb[ i ] < this->get_dbl_par(dblRelAcc)) {
	    EmptySet = true;
	    break;
	    }

	   if( Lpt[ i ] > ub[ i ] ) {
	    Lpt[ i ] = ub[ i ];
	    LptHasChgd = true;
	    break;
	    }
	   else
	    if( Lpt[ i ] < lb[ i ] ) {
	     Lpt[ i ] = lb[ i ];
	     LptHasChgd = true;
	     break;
	     }
	   }

	  if( ! EmptySet )
	   for( ; i < NumVar ; i++ ) {
	    if( ub[ i ] - lb[ i ] < this->get_dbl_par(dblRelAcc) ) {
	     EmptySet = true;
	     break;
	     }

	    if( Lpt[ i ] > ub[ i ] )
	     Lpt[ i ] = ub[ i ];
	    else
	     if( Lpt[ i ] < lb[ i ] )
	      Lpt[ i ] = lb[ i ];
	    }
	  }
	 else
	  if( lb.empty() ) {
	   Variable::Index i = 0;
	   for( ; i < NumVar ; i++ )
	    if( Lpt[ i ] < lb[ i ] ) {
	     Lpt[ i ] = lb[ i ];
	     LptHasChgd = true;
	     }
	   for( ; i < NumVar ; i++ )
	    if( Lpt[ i ] < lb[ i ] )
	     Lpt[ i ] = lb[ i ];
	   }
	  else {
		  Variable::Index i = 0;
	   for(  ; i < NumVar ; i++ )
	    if( Lpt[ i ] > ub[ i ] ) {
	     Lpt[ i ] = ub[ i ];
	     LptHasChgd = true;
	     }
	   for( ; i < NumVar ; i++ )
	    if( Lpt[ i ] > ub[ i ] )
	     Lpt[ i ] = ub[ i ];
	   }

	 return( LptHasChgd );
}

void SubGradSolver::ProjTng(std::vector<double> Gpt, std::vector<double> Lpt,
		const Variable::Index strt, Variable::Index stp) {

	// project Gpt in the tangent cone at Lpt  - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( stp > NumVar )
	  stp = NumVar;

	 #if SubGrad_HANDLES_CONSTRAINTS
	 #endif

	 if( ( ! lb.empty() ) && ( ! ub.empty() ) )
		 return;

	 // project Gpt on the box constraints  - - - - - - - - - - - - - - - - - -

	  if( lb.empty() && ub.empty() ) {
	   for( Variable::Index j = strt ; j < stp ; j++ )
	    if( ( ( Lpt[ j ] >= ub[ j ] ) && ( Gpt[ j ] <= 0 ) ) ||
	        ( ( Lpt[ j ] <= lb[ j ] ) && ( Gpt[ j ] >= 0 ) ) )
	     Gpt[ j ] = 0;
	   }
	  else
	   if( lb.empty() ) {
	    for( Variable::Index j = strt ; j < stp ; j++ )
	     if( ( Lpt[ j ] <= lb[ j ] ) && ( Gpt[ j ] >= 0 ) )
	      Gpt[ j ] = 0;
	    }
	   else
	    for( Variable::Index j = strt ; j < stp ; j++ )
	     if( ( Lpt[ j ] >= ub[ j ] ) && ( Gpt[ j ] <= 0 ) )
	      Gpt[ j ] = 0;
}

void SubGradSolver::InitStepsize(void) {
}

void SubGradSolver::InitDeflection(void) {
}

void SubGradSolver::ChgLambdaHat(void) {

	 // compute the LambdaHat- - - - - - - - - - - - - - - -  - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( alpha == 1 ) {
	  HatEpsilon = SigmaHat;
	  LambdaHat.swap(Lambda);
	  FiHat = FiLambda;
	  }
	 else {
	  // compute the linearization error of DirM1 in LambdaHat- - - - - - - - - -
	  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  double SclPrd = 0;
	  for( Variable::Index i = 0 ; i < NumVar ; i++ )
	   SclPrd += dir[ i ] * ( LambdaHat[ i ] - Lambda[ i ] );

	  HatEpsilon += alpha * ( FiLambda - FiHat + SclPrd );

	  // compute the linearization error of Gi[wFi] in LambdaHat  - - - - - - - -
	  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  SclPrd = 0;
	  for( Variable::Index i = 0 ; i < NumVar ; i++ )
	   SclPrd += Gi[ i ] * ( Lambda[ i ] - LambdaHat[ i ] );

	  SigmaHat += (1.0 - alpha ) * ( FiHat - FiLambda + SclPrd );

	  // compute the LambdaHat and FiHat - - - - - - - - - -  - - - - - - - - - -
	  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  // g := x1 * g1 + x2 * g2
	  Function::Vec_FunctionValue p_v,p_v1;
	  p_v.swap(Lambda);
	  p_v1.swap(LambdaHat);

	  std::transform(p_v.begin(), p_v.end(), p_v.begin(),
					 std::bind1st(std::multiplies<double>(), alpha)); // x1 * g1
	  std::transform(p_v1.begin(), p_v1.end(), p_v1.begin(),
							 std::bind1st(std::multiplies<double>(), double(1)-alpha)); //x2 * g2
	  std::transform (p_v.begin(), p_v.end(), p_v1.begin(), LambdaHat.begin(), std::plus<double>());



	  FiHat = alpha * FiLambda + (1.0 - alpha) * FiHat;

	  // compute the linearization error of Dir in LambdaHat  - - - - - - - - - -
	  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  HatEpsilon = alpha * SigmaHat + ( double(1) - alpha ) * HatEpsilon;
	  }  // end( alpha < 1 )- - - - - - - - - - - - - - - - - - - - - - - - - - -

	 VLOG( 1 , " ~ HFi = " << -FiHat );



}

void SubGradSolver::UpdateSigma(void) {

	// some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( ! deflection )
	  throw( std::logic_error( "SubGrad::UpdateSigma: this is not allowed" ) );

	 // compute the linearization error of Gi[wFi] at Lambdabar - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 Sigma += FiBar - FiLambda;
	 if( ! LHasProj )
	  Sigma -= step * std::inner_product(Gi.begin(),Gi.end(),dir.begin(),0);// ScalarProduct( Gi , dir , NumVar );
	 else {
	  double SclPrd = 0;
	  for( Variable::Index i = 0 ; i < NumVar ; i++ )
	   SclPrd +=  Gi[ i ] * ( Lambda[ i ] - LambdaBar[ i ] );
	  Sigma += SclPrd;
	  }

	 VLOG( 1 , " ~ Sigma1 = " << Sigma );



}

//void SubGradSolver::UpdtLowerBound(void) {
//	 const double  LwrBnd = Oracle->GetLowerBound();
//	 if( LwrBnd > LowerBound ) {
//	  LowerBound = LwrBnd;
//	  TrueLB = ( LwrBnd > Oracle->GetMinusInfinity() );
//	  }
//}

void SubGradSolver::Log1(void) {

//#if( SubGrad_LOG )
  if( LogVerb > 1 ) {
   *f_log << std::endl << "{" << SCalls << "-" << ParIter << "} ";

   *f_log << " FiC = ";
    if( FiBar == Inf<double>() )
     *f_log << " -INF";
    else
      *f_log << -FiBar;

    *f_log << " ~ FiBest = ";
    if( FiBest == Inf<double>() )
     *f_log << " -INF";
    else {
     *f_log << -FiBest;

     if( LowerBound > - Inf<double>() )
      *f_log << " ~ Gap = " << ( FiBest - LowerBound ) /
                                std::max( std::abs( LowerBound ) , 1.0 ) << " %";
     }
   }
}

void SubGradSolver::Log2(void) {
//#if( SubGrad_LOG )
 // information relative to the direction

 if( LogVerb > 1 ) {
  *f_log << std::endl << "           ";

  if( KSGPar4 )
   *f_log << " ( |D| , Epsilon , HatEpsilon ) = (" << sqrt(NrmDir)
	    << " , " << Epsilon << " , " << HatEpsilon <<  ") ~ GixD = "
	    << dGk;
  else
   *f_log << " ( |D| , Epsilon ) = (" << sqrt(NrmDir) << " , "
	    << Epsilon <<  ") ~ GixD = " << dGk;

  }
//#endif
}

void SubGradSolver::Log3(void) {
//#if( SubGrad_LOG )
  if( LogVerb > 1 ) {
   *f_log << std::endl << "           ";

   *f_log << " Fi1 = ";  // print function value
   if( FiLambda == Inf<double>() )
     *f_log << " -INF => STOP.";
   else
    *f_log << -FiLambda;

   }
 //#endif
}

void SubGradSolver::Log4(void) {
//#if( SubGrad_LOG )
  if( LogVerb > 1 ) {
   // if the deflection is used the information relative to the subgradient
   // have already been displayed [see Log2() ]

   *f_log << std::endl << "           ";
//   if( deflection)
//    *f_log << " ( |Gi| , Sigma ) = (" << sqrt(NrmGi) << " , " << Sigma
//	    << ")";
   }
// #endif
}

void SubGradSolver::FormD(void) {

	// project the subgradient and compute its norm - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 // take in account some comments and details about FormD():
	 // (i) the method computes the sugradient in Lambda, if not already computed;
	 // (ii) the method computes the deflection coefficient Alpha, and a fortiori
	 //      the direction Dir;
	 // (iii) the method computes the stepsize Nu using the deflection-restricted
	 //       scheme;
	 // (iv) if the variables-set has been changed after a call to GetFiStatus,
	 //      even the projection of both Gi and Dir is done in AddVariables or
	 //      RemoveVariables.
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( ( fs == kFiNorm ) && ( ! DoSS ) && ( KSGPar1 & PrjSG ) ) {
	  // the projection of Gi[wFi] is already performed but relatively to the
	  // tangent cone at the tentative point of the previous iteration: if a NS
	  // occurred, the projection has to be recomputed


	  //Variable::Index SGBDim = Oracle->GetGi( Gi , SGBse );
	  //Variable::Index SGBDim;
	  //const C05Function::LinearizationName SGBDim;

	 // cIndex_Set SGBse;
	   // Index SGBDim = Oracle->GetGi( Gi , SGBse );


	  std::vector<double> pg(NumVar);
	  for(auto c05 : v_c05f){
		  c05->get_linearization_coefficients(pg.data());
		  std::transform (Gi.begin(), Gi.end(), pg.begin(), Gi.begin(), std::plus<double>());


	  }
	  // if the case of when the 0-th component is a linear function
	  if(lf){
		  lf->get_linearization_coefficients(pg.data());
		  std::transform (Gi.begin(), Gi.end(), pg.begin(), Gi.begin(), std::plus<double>());

	  }




	  ProjTng( Gi , LambdaBar );
	  NrmGi = std::inner_product(Gi.begin(), Gi.end(), Gi.begin(), 0);


	  dM1GkDone = false;
	  }

	 // subgradient log: print the information relative to Gi[wFi] - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 Log4();

	 // project DirM1 and compute its norm - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if ( ( KSGPar1 & PrjDIRM1 ) && ( fs == kFiNorm ) )  {
	  if( KSGPar1 & PrjDIR ) {  // DirM1 is projected over the tangent cone at the
		                   // previous LambdaBar - - - - - - - - - - - - - - -
	   if ( DoSS )  {  // the projection is no longer valid, re-project DirM1
	    ProjTng( dir , LambdaBar );
	    NrmDir = std::inner_product(dir.begin(), dir.end(), dir.begin(), 0);
	    dM1GkDone = false;
	    }
	   }
	  else {  // project DirM1 - - - - - - - - - - - - - - - - - - - - - - - - - -
	   ProjTng( dir , LambdaBar );
	   NrmDir = std::inner_product(dir.begin(), dir.end(), dir.begin(), 0);

	   dM1GkDone = false;
	   }
	  }  // end DirM1 projection - - - - - - - - - - - - - - - - - - - - - - - - -
	     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 // (re)compute the scalar product between DirM1 and Gi[wFi] at Lambda - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( deflection) {
	  if( ! dM1GkDone )
		dM1Gk = std::inner_product(Gi.begin(),Gi.end(), dir.begin(),0);

	  }
	 else {  // if no deflection is done, Dir <==> Gi[wFi]
	  dGk = NrmDir = NrmGi;
	  if( ! dM1GkDone )
	   dM1Gk = dGk;
	  }

	 dM1GkDone = true;

	 // stepsise and deflection computation  - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( KSGPar3 & StepRST ) {  // stepsize-restricted approach - - - - - - - - - -

	  if( deflection) {  // deflection coefficient computation
	   deflection->NewDEF(  );
	   alpha = deflection->GetDFLCoeff();
	   }

	  VLOG( 1 , " ~ alpha = " << alpha );
	  }
	 else {  // deflection-restricted approach - - - - - - - - - - - - - - - - - -
	         //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  dGk = dM1Gk;  // this scheme doesn't distinguish between dGk and dM1Gk

	  if( ! InnIter )  // stepsize computation - - - - - - - - - - - - - - - - - -
	   stepsize->NewStep( );

	  step = stepsize->GetStepsize( InnIter );
	  VLOG( 1 , " ~ step = " << step );

	  deflection->NewDEF( );  // no deflection is not allowed
	  alpha = deflection->GetDFLCoeff();
	  VLOG( 1 , " ~ alpha = " << alpha );

	  // safe rule - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  if( ( KSGPar3 & SafeRULE )  && !InnIter ) {
	   double target =  ( step * NrmDir ) /
	                   ( FiBar - stepsize->GetLev() + step * NrmDir );
	   alpha = std::min( 1.0 , std::max( alpha , target ) );
	   VLOG( 1 , " ~ fixed alpha = " << alpha );
	   }
	  }  // end (deflection-restricted approach) - - - - - - - - - - - - - - - - -
	     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 // tell the oracle to keep track the direction- - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 // DirPos == first => DirM1 is at the position MaxNCnst + 1
	 // and Dir is registered in the free entry MaxNCnst + 2,
	 // and vice versa for the case DirPos == second


	 if( deflection ) {
	  if(  alpha == 1 ) {  // pass the subgradient to the FiOracle object

	   C05Function::LinearCombination lc({std::make_pair(0,alpha)});




	   for(auto c05 : v_c05f)
		   c05->store_combination_of_linearizations(lc , DirPos == first ? 2 : 1);

	   if(lf)
		   lf->store_combination_of_linearizations(lc , DirPos == first ? 2 : 1);
	   }
	  else {

	   C05Function::LinearCombination lc(2);
	   lc={{1,alpha},{2,1.0 - alpha }};

	   if( DirPos == first ) {

	    for(auto c05 : v_c05f){
	    	   		  c05->store_combination_of_linearizations(lc, 2);
	    	   }

	    if(lf){
	    	   		  lf->store_combination_of_linearizations(lc , 2);


	    	   	  }

	    }
	   else {
		   for(auto c05 : v_c05f){
		   	   		  c05->store_combination_of_linearizations(lc, 1);
		   	   }

		   if(lf){
		   	   	lf->store_combination_of_linearizations(lc , 1);

		    }

	    }
	   }
	  }

	 // print Dir norm - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 VLOG( 1 , " ~ trial |D| = " << NrmDir );


}

//void SubGradSolver::MemDealloc(void) {
//
//}

void SubGradSolver::FormLambda1(void) {

	if( deflection )
	  for( Variable::Index i = 0 ; i < NumVar ; i++ )
	   Lambda[ i ] = LambdaBar[ i ] - step * dir[ i ];
	 else
	  for( Variable::Index i = 0 ; i < NumVar ; i++ )
	   Lambda[ i ] = LambdaBar[ i ] - step * Gi[ i ];

	 LHasChgd = true;

}

Solver::OFValue SubGradSolver::get_lb(void) {
}

Solver::OFValue SubGradSolver::get_ub(void) {
}

void SubGradSolver::FiAndGi(const Variable::Index wFi) {

	// pass Lambda to the Oracle - - - - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( LHasChgd ) {  // if Lambda1 has changed, set the new one to the Oracle


	  int idl=-1;
	  for(auto l : LamVcblr){
		  auto lc= dynamic_cast<ColVariable*>(l);
		  lc->set_value(Lambda[idl++]);
	  }

	  //Oracle->SetLambda( Lambda );
	  LHasChgd = false;
	  }

	 // if outer iteration, evaluate Fi in Lambda1  - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( wFi == Inf<Variable::Index>() ) {  // call Fi()- - - - - - - - - - - - - - - - - -

	  for(int i=0; i<v_c05f.size();i++){
		  v_c05f[i]->compute(true);
		  FiLambda+=v_c05f[i]->get_value();

	  }
	  if(lf){
		  FiLambda+=lf->get_value();
	  }

	  //FiLambda = Oracle->Fi( );
	  FiEvaltns++;

	  if( FiLambda == -Inf<double>() )    // Fi() unbounded
	   return;

	  // update FiBest, if necessary - - - - - - - - - - - - - - - - - - - - - - -
	  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  if( FiLambda < FiBest ) {
	   FiBest = FiLambda;
	   if( KpBstL )
		   LambdaBest.swap(Lambda);

	   }
	  }  // end( function evaluation - - - - - - - - - - - - - - - - - - - - - - -

	 // if Lambda is unfeasible, a new subgradient is computed; otherwise a new
	 // constraint is added - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 //if( FiLambda == Inf<double>() )  // it is a constraint
	 // #if SubGrad_HANDLES_CONSTRAINTS
	 // #endif
	 //else
	  EvalGi( wFi );  // the item is a subgradient



}



void SubGradSolver::GotoLambda1(void) {

	// compute the linearization error of Dir at the new LambdaBar  - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( FiBar < Inf<double>() && deflection ) {
	  double SclPrd;
	  if( ! LHasProj )
	   SclPrd = step * NrmDir;
	  else {
	   SclPrd = 0;
	   for( Variable::Index i = 0 ; i < NumVar ; i++ )
	    SclPrd += dir[ i ] * ( LambdaBar[ i ] - Lambda[ i ] );
	   }

	  Epsilon += FiLambda - FiBar + SclPrd;
	  }

	 // move Lambda to LambdaBar= - - - - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 LambdaBar.swap(Lambda);
	 //VectAssign( LambdaBar , Lambda , NumVar );
	 if( ! InnIter )
	  FiBar = FiLambda;


}

void SubGradSolver::MemDealloc(void) {

	// the point - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 delete[] LambdaBar;
	 delete[] Lambda;
	 delete[] LambdaBest;
	 delete[] LambdaHat;

	 // box and simplex constraints - - - - - - - - - - - - - - - - - - - - - - -

	 delete[] lb;
	 delete[] ub;

	 #if SubGrad_HANDLES_CONSTRAINTS
	 #endif

	 // the subgradients and direction- - - - - - - - - - - - - - - - - - - - - -

	 delete[] Gi;
	 delete[] dir;

	 delete[] MultBse;

}

void SubGradSolver::SaveDir(void) {

	// compute the linearization error of Dir at LambdaBar- - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( alpha == 1 )
	  Epsilon = Sigma;
	 else
	  Epsilon = alpha * Sigma + ( double( 1 ) - alpha ) * Epsilon;

	 // change LambdaHat - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( KSGPar4 )
	  ChgLambdaHat();

	 // save the search direction Dir  - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( deflection ) {
	  // Dir is the convex combination of DirM1 and Gi[wFi]- - - - - - - - - - - -
	  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  if( alpha == 1 )
		  dir.swap(Gi);

	  else{
		  // g := x1 * g1 + x2 * g2
		  Function::Vec_FunctionValue p_v,p_v1;
		  p_v.swap(Gi);
		  p_v1.swap(dir);

		  std::transform(p_v.begin(), p_v.end(), p_v.begin(),
						 std::bind1st(std::multiplies<double>(), alpha)); // x1 * g1
		  std::transform(p_v1.begin(), p_v1.end(), p_v1.begin(),
								 std::bind1st(std::multiplies<double>(), double(1)-alpha)); //x2 * g2
		  std::transform (p_v.begin(), p_v.end(), p_v1.begin(), dir.begin(), std::plus<double>());


	  }

	  // swap the position for DirM1 and Dir - - - - - - - - - - - - - - - - - - -
	  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  if( DirPos == first )
	   DirPos = second;
	  else
	   DirPos = first;

	  // project Dir on tangent cone at LambdaBar  - - - - - - - - - - - - - - - -
	  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  if( KSGPar1 & PrjDIR )
	   ProjTng( dir , LambdaBar );

	  // compute the (squared) norm of the direction and its scalar product with
	  // the subgradient   - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  NrmDir = std::inner_product(dir.begin(),dir.end(),dir.begin(),0);

	  dGk = std::inner_product(dir.begin(),dir.end(),Gi.begin(),0);


	  }  // end( deflection )  - - - - - - - - - - - - - - - - - - - - - - - - - -
	     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 // stepsize computation - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( KSGPar3 & StepRST ) {
	  if( ! InnIter )
	   stepsize->NewStep( );

	  step = stepsize->GetStepsize( InnIter );
	  VLOG( 1 , " ~ step = " << step );

	  // safe rule - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  if( KSGPar3 & SafeRULE ) {
	   stepsize->SetMaxBeta( alpha );
	   step = stepsize->GetStepsize( InnIter );
	   VLOG( 1 , " ~ (fixed) step = " << step );
	   }

	  }  // end stepsize-restricted scheme - - - - - - - - - - - - - - - - - - - -


}
