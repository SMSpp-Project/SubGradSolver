/*--------------------------------------------------------------------------*/
/*------------------------ File SubGradSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SubGradSolver class, which implements the CDASolver
 * interface by means of a (deflected, projected, incremental) subgradient
 * method, and of the methods of the base classes Stepsize and Deflection
 * that give the rules access to the data of SubGradSolver.
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SubGradSolver.h"

#include "Block.h"

#include "DQuadFunction.h"

#include "FRealObjective.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include <cmath>

#include <iomanip>

#include <numeric>

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define SGLOG( l , x ) if( f_log && ( LogVerb > l ) ) *f_log << x

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using Index = SubGradSolver::Index;
using VarValue = SubGradSolver::VarValue;
using Vec_VarValue = SubGradSolver::Vec_VarValue;
using LinearCombination = SubGradSolver::LinearCombination;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr int PrjSG = 1;     // project g_i
static constexpr int PrjDIRM1 = 2;  // project d_{i-1}
static constexpr int PrjDIR = 4;    // project d_i

static constexpr int SafeRULE = 1;  // the safe rule
static constexpr int StepRST = 2;   // the stepsize-restricted scheme

static constexpr int RstRules = 1;   // reset the rules at each compute()
static constexpr int RstCenter = 2;  // re-read the center at each compute()

// each component has kPoolSz names in its global pool, appended to those
// that are there when SubGradSolver is attached, so as not to touch those
// of other Solver: the latest linearization is in the first, the direction
// alternately in the other two, the new one being formed out of the old one
static constexpr Index kSlotG = 0;
static constexpr Index kPoolSz = 3;

// the method stalls if the stepsize stays below SmallStep * dbltStar for
// more than MaxSmallStep consecutive (outer) iterations
static constexpr double SmallStep = 1e-8;
static constexpr Index MaxSmallStep = 100;

static constexpr VarValue INF = Inf< VarValue >();

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( SubGradSolver );

/*--------------------------------------------------------------------------*/
/*------------------------- PARAMETER NAME TABLES --------------------------*/
/*--------------------------------------------------------------------------*/

static const std::vector< std::string > SGS_int_pars_str = {
 "intProjection" , "intScheme" , "intKeepHat" , "intSeed" , "intResetAlg" ,
 "intPrimalRec" , "intColorGreen" , "intColorYellow" , "intColorRed" ,
 "intFumeroR1" , "intFumeroEta1" , "intFumeroEta2" , "intVolumeTauIter" ,
 "intPDAverage"
 };

static const std::vector< std::string > SGS_dbl_pars_str = {
 "dbltStar" , "dblNZEps" , "dblIncrFactor" , "dblLipFactor" , "dblPolyakBeta" ,
 "dblColorBeta0" , "dblFumeroBeta0" , "dblFumeroSigmaMin" ,
 "dblVolumeTauInit" , "dblVolumeTauMin" , "dblVolumeTauFactor" ,
 "dblVolumeM"
 };

static const std::vector< std::string > SGS_str_pars_str = {
 "strStepsize" , "strDeflection"
 };

static const std::vector< int > SGS_dflt_int_par = {
 0 ,    // intProjection
 1 ,    // intScheme
 1 ,    // intKeepHat
 0 ,    // intSeed
 0 ,    // intResetAlg
 0 ,    // intPrimalRec
 1 ,    // intColorGreen
 400 ,  // intColorYellow
 10 ,   // intColorRed
 5 ,    // intFumeroR1
 5 ,    // intFumeroEta1
 5 ,    // intFumeroEta2
 100 ,  // intVolumeTauIter
 1      // intPDAverage
 };

static const std::vector< double > SGS_dflt_dbl_par = {
 100 ,   // dbltStar
 0 ,     // dblNZEps
 0 ,     // dblIncrFactor
 1 ,     // dblLipFactor
 1e-2 ,  // dblPolyakBeta
 1 ,     // dblColorBeta0
 1 ,     // dblFumeroBeta0
 1e-4 ,  // dblFumeroSigmaMin
 1 ,     // dblVolumeTauInit
 1e-4 ,  // dblVolumeTauMin
 0.5 ,   // dblVolumeTauFactor
 0.1     // dblVolumeM
 };

static const std::vector< std::string > SGS_dflt_str_par = {
 "ColorTV" ,  // strStepsize
 "Volume"     // strDeflection
 };

// the names, the old ones included
static const std::map< std::string , ThinComputeInterface::idx_type >
SGS_int_pars_map = {
 { "intProjection" , SubGradSolver::intProjection } ,
 { "intScheme" , SubGradSolver::intScheme } ,
 { "intKeepHat" , SubGradSolver::intKeepHat } ,
 { "intSeed" , SubGradSolver::intSeed } ,
 { "intResetAlg" , SubGradSolver::intResetAlg } ,
 { "intPrimalRec" , SubGradSolver::intPrimalRec } ,
 { "intColorGreen" , SubGradSolver::intColorGreen } ,
 { "intColorYellow" , SubGradSolver::intColorYellow } ,
 { "intColorRed" , SubGradSolver::intColorRed } ,
 { "intFumeroR1" , SubGradSolver::intFumeroR1 } ,
 { "intFumeroEta1" , SubGradSolver::intFumeroEta1 } ,
 { "intFumeroEta2" , SubGradSolver::intFumeroEta2 } ,
 { "intVolumeTauIter" , SubGradSolver::intVolumeTauIter } ,
 { "intPDAverage" , SubGradSolver::intPDAverage }
 };

static const std::map< std::string , ThinComputeInterface::idx_type >
SGS_dbl_pars_map = {
 { "dbltStar" , SubGradSolver::dbltStar } ,
 { "dblNZEps" , SubGradSolver::dblNZEps } ,
 { "dblIncrFactor" , SubGradSolver::dblIncrFactor } ,
 { "dblLipFactor" , SubGradSolver::dblLipFactor } ,
 { "dblPolyakBeta" , SubGradSolver::dblPolyakBeta } ,
 { "dblColorBeta0" , SubGradSolver::dblColorBeta0 } ,
 { "dblFumeroBeta0" , SubGradSolver::dblFumeroBeta0 } ,
 { "dblFumeroSigmaMin" , SubGradSolver::dblFumeroSigmaMin } ,
 { "dblVolumeTauInit" , SubGradSolver::dblVolumeTauInit } ,
 { "dblVolumeTauMin" , SubGradSolver::dblVolumeTauMin } ,
 { "dblVolumeTauFactor" , SubGradSolver::dblVolumeTauFactor } ,
 { "dblVolumeM" , SubGradSolver::dblVolumeM }
 };

static const std::map< std::string , ThinComputeInterface::idx_type >
SGS_str_pars_map = {
 { "strStepsize" , SubGradSolver::strStepsize } ,
 { "strDeflection" , SubGradSolver::strDeflection }
 };

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

static inline VarValue norm2( const Vec_VarValue & v )
{
 return( std::inner_product( v.begin() , v.end() , v.begin() ,
			     VarValue( 0 ) ) );
 }

/*--------------------------------------------------------------------------*/

static inline VarValue dot( const Vec_VarValue & v , const Vec_VarValue & w )
{
 return( std::inner_product( v.begin() , v.end() , w.begin() ,
			     VarValue( 0 ) ) );
 }

/*--------------------------------------------------------------------------*/
// v = a * v + b * w

static inline void comb( Vec_VarValue & v , VarValue a ,
			 const Vec_VarValue & w , VarValue b )
{
 for( Index i = 0 ; i < v.size() ; ++i )
  v[ i ] = a * v[ i ] + b * w[ i ];
 }

/*--------------------------------------------------------------------------*/
// true if s names "no Deflection"

static inline bool is_none( const std::string & s )
{
 return( s.empty() || ( s == "none" ) || ( s == "None" ) );
 }

/*--------------------------------------------------------------------------*/
// scientific notation with the given number of digits

static inline std::ostream & sci( std::ostream & os , int prec )
{
 os.setf( std::ios::scientific , std::ios::floatfield );
 os << std::setprecision( prec );
 return( os );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF SubGradSolver -------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------- CONSTRUCTOR AND DESTRUCTOR ------------------------*/
/*--------------------------------------------------------------------------*/

SubGradSolver::SubGradSolver( void )
 : CDASolver() , f_convex( true ) , f_zero( nullptr ) , NrFi( 0 ) ,
   NumVar( 0 ) , f_empty_box( false ) , stepsize( nullptr ) ,
   deflection( nullptr ) , f_shared_rule( false ) , f_rules_dirty( true ) ,
   f_scheme( 0 ) , Result( kUnEval ) , SCalls( 0 ) , ParIter( 0 ) ,
   TotIter( 0 ) , FiEvaltns( 0 ) , GiEvaltns( 0 ) , FiBar( INF ) ,
   FiLambda( INF ) , FiBest( INF ) , FiHat( INF ) , GiCnst( 0 ) ,
   dirCnst( 0 ) , alpha( 1 ) , step( 0 ) , Sigma( 0 ) , Epsilon( 0 ) ,
   SigmaHat( 0 ) , HatEpsilon( 0 ) , NrmGi( 0 ) , NrmDir( 0 ) ,
   NrmPrj( INF ) , dGk( 0 ) , dM1Gk( 0 ) , dM1GkDone( false ) ,
   MaxNrmGi( 0 ) , LHasProj( false ) ,
   DoSS( true ) , f_dir_valid( false ) , f_eval_ok( false ) , DirPos( 1 ) ,
   TrueLB( - INF ) , CondLB( - INF ) , CertLB( - INF ) , CSSCntr( 0 ) ,
   CNSCntr( 0 ) , CSmallStep( 0 ) , NItIncr( 0 ) , InnIter( false ) ,
   f_need_rebuild( false )
{
 MaxIter = get_dflt_int_par( intMaxIter );
 LogVerb = get_dflt_int_par( intLogVerb );
 EverykIt = get_dflt_int_par( intEverykIt );
 SGPar1 = get_dflt_int_par( intProjection );
 SGPar3 = get_dflt_int_par( intScheme );
 SGPar4 = get_dflt_int_par( intKeepHat );
 SGPar5 = get_dflt_int_par( intSeed );
 RstAlg = get_dflt_int_par( intResetAlg );
 PrimalRec = get_dflt_int_par( intPrimalRec );

 MaxTime = get_dflt_dbl_par( dblMaxTime );
 EveryTTm = get_dflt_dbl_par( dblEveryTTm );
 RelAcc = get_dflt_dbl_par( dblRelAcc );
 AbsAcc = get_dflt_dbl_par( dblAbsAcc );
 UpCutOff = get_dflt_dbl_par( dblUpCutOff );
 LwCutOff = get_dflt_dbl_par( dblLwCutOff );
 tStar = get_dflt_dbl_par( dbltStar );
 NZEps = get_dflt_dbl_par( dblNZEps );
 SGPar2 = get_dflt_dbl_par( dblIncrFactor );

 StepsizeName = get_dflt_str_par( strStepsize );
 DeflectionName = get_dflt_str_par( strDeflection );

 for( idx_type p = intColorGreen ; p < intLastParSGS ; ++p )
  v_rule_int.push_back( get_dflt_int_par( p ) );

 for( idx_type p = dblLipFactor ; p < dblLastParSGS ; ++p )
  v_rule_dbl.push_back( get_dflt_dbl_par( p ) );

 ErgSum = 0;
 myrandom.seed( SGPar5 );

 c_start = std::chrono::system_clock::now();

 v_events.resize( max_event_number() );
 }

/*--------------------------------------------------------------------------*/

SubGradSolver::~SubGradSolver()
{
 set_Block( nullptr );
 delete_rules();
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void SubGradSolver::set_Block( Block * block )
{
 if( f_Block == block )  // registering to the same Block
  return;                // cowardly and silently return

 if( f_Block )           // changing from a previous Block
  clear_structure();

 Solver::set_Block( block );  // attach to the new Block

 if( ! f_Block )  // that was actually clearing the Block
  return;         // all done

 // the abstract representation may have to be generated, hence the Block is
 // lock()-ed rather than read_lock()-ed
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
		     "SubGradSolver::set_Block: unable to lock the Block" ) );

 try {
  build_structure( true );
  }
 catch( ... ) {
  if( ! owned )
   f_Block->unlock( f_id );
  throw;
  }

 if( ! owned )
  f_Block->unlock( f_id );

 // the rules are constructed now if their parameters have been set already,
 // as when a BlockSolverConfig sets them before registering the Solver
 if( f_rules_dirty || ( ! stepsize ) )
  build_rules();

 reset_function_info();

 }  // end( SubGradSolver::set_Block )

/*--------------------------------------------------------------------------*/

void SubGradSolver::build_structure( bool generate )
{
 if( generate ) {
  f_Block->generate_abstract_variables();
  f_Block->generate_abstract_constraints();
  f_Block->generate_objective();
  }

 // the sub-Block that this Solver has been told to ignore are not
 // components [see Solver::set_excluded_blocks()]
 std::vector< Block * > sb;
 for( auto blk : f_Block->get_nested_Blocks() )
  if( ! is_excluded( blk ) )
   sb.push_back( blk );

 v_c05f.clear();
 f_zero = nullptr;

 if( sb.empty() ) {  // no sub-Block: the objective is a C05Function - - - - -
  auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
  if( ! obj )
   throw( std::invalid_argument( "SubGradSolver::set_Block: "
				 "objective is not a FRealObjective" ) );

  auto c05f = dynamic_cast< C05Function * >( obj->get_function() );
  if( ! c05f )
   throw( std::invalid_argument( "SubGradSolver::set_Block: "
				 "the objective is not a C05Function" ) );

  f_convex = c05f->is_convex();
  if( ( ! f_convex ) && ( ! c05f->is_concave() ) )
   throw( std::invalid_argument( "SubGradSolver::set_Block: "
				 "only convex or concave objectives allowed"
				 ) );

  if( ( f_convex && ( obj->get_sense() == Objective::eMax ) ) ||
      ( ( ! f_convex ) && ( obj->get_sense() == Objective::eMin ) ) ) {
   // a LinearFunction is both convex and concave: what counts is the sense
   if( c05f->is_convex() && c05f->is_concave() )
    f_convex = ( obj->get_sense() == Objective::eMin );
   else
    throw( std::invalid_argument( "SubGradSolver::set_Block: "
				  "can only minimize convex / maximize "
				  "concave" ) );
   }

  v_c05f.push_back( c05f );
  }
 else {  // there are sub-Block: each has a C05Function objective- - - - - -
  if( auto obj = f_Block->get_objective() ) {
   auto fobj = dynamic_cast< FRealObjective * >( obj );
   if( ! fobj )
    throw( std::logic_error( "SubGradSolver::set_Block: "
			     "the objective is not a real function" ) );

   if( auto fun = fobj->get_function() ) {
    if( ( ! dynamic_cast< LinearFunction * >( fun ) ) &&
	( ! dynamic_cast< DQuadFunction * >( fun ) ) )
     throw( std::logic_error( "SubGradSolver::set_Block: the objective is "
			      "neither a LinearFunction nor a "
			      "DQuadFunction" ) );
    f_zero = dynamic_cast< C05Function * >( fun );
    if( ! f_zero->get_num_active_var() )
     f_zero = nullptr;  // no Variable: no 0-th component
    }
   }

  for( Index i = 0 ; i < sb.size() ; ++i ) {
   auto obj = dynamic_cast< FRealObjective * >( sb[ i ]->get_objective() );
   if( ! obj )
    throw( std::logic_error( "SubGradSolver::set_Block: sub-Block "
			     "objective is not a real function" ) );

   auto c05f = dynamic_cast< C05Function * >( obj->get_function() );
   if( ! c05f )
    throw( std::logic_error( "SubGradSolver::set_Block: sub-Block "
			     "objective is not a C05Function" ) );

   if( ! i )
    f_convex = c05f->is_convex();
   else
    if( c05f->is_convex() != f_convex )
     throw( std::invalid_argument( "SubGradSolver::set_Block: objectives "
				   "must be all convex or all concave" ) );

   if( ( ! f_convex ) && ( ! c05f->is_concave() ) )
    throw( std::invalid_argument( "SubGradSolver::set_Block: only convex "
				  "or concave objectives allowed" ) );

   if( ( f_convex && ( obj->get_sense() == Objective::eMax ) ) ||
       ( ( ! f_convex ) && ( obj->get_sense() == Objective::eMin ) ) )
    throw( std::invalid_argument( "SubGradSolver::set_Block: can only "
				  "minimize convex / maximize concave" ) );

   if( ! sb[ i ]->get_nested_Blocks().empty() )
    throw( std::logic_error( "SubGradSolver::set_Block: "
			     "nephew sub-Blocks are not allowed" ) );

   for( auto groups : { & sb[ i ]->get_static_variable_groups() ,
			& sb[ i ]->get_dynamic_variable_groups() ,
			& sb[ i ]->get_static_constraint_groups() ,
			& sb[ i ]->get_dynamic_constraint_groups() } )
    for( const auto & group : *groups )
     if( group && group->get_num_elements() )
      throw( std::logic_error( "SubGradSolver::set_Block: Variable or "
			       "Constraint in sub-Block are not allowed" ) );

   v_c05f.push_back( c05f );
   }
  }

 NrFi = v_c05f.size();

 // the ColVariable: the union of the active Variable of the functions, the
 // 0-th first, in the order they are first met; the position in LamVcblr of
 // each active Variable of each function is in v_loc2glob, which is left
 // empty when it is the identity - - - - - - - - - - - - - - - - - - - - - -

 LamVcblr.clear();
 v_loc2glob.assign( NrFi + 1 , {} );
 std::unordered_map< ColVariable * , Index > var2idx;

 auto scan = [ & ]( C05Function * f , std::vector< Index > & map ) {
  map.reserve( f->get_num_active_var() );
  for( auto it = std::as_const( f )->begin() ;
       it != std::as_const( f )->end() ; ++it ) {
   auto p = dynamic_cast< ColVariable * >( & ( *it ) );
   if( ! p )
    throw( std::logic_error( "SubGradSolver::set_Block: some active "
			     "Variable is not a ColVariable" ) );
   auto [ jt , ins ] = var2idx.try_emplace( p , LamVcblr.size() );
   if( ins )
    LamVcblr.push_back( p );
   map.push_back( jt->second );
   }
  };

 if( f_zero )
  scan( f_zero , v_loc2glob[ NrFi ] );

 for( Index k = 0 ; k < NrFi ; ++k )
  scan( v_c05f[ k ] , v_loc2glob[ k ] );

 NumVar = LamVcblr.size();

 for( auto & map : v_loc2glob ) {
  bool id = ( map.size() == NumVar );
  for( Index i = 0 ; id && ( i < map.size() ) ; ++i )
   if( map[ i ] != i )
    id = false;
  if( id )
   map.clear();
  else
   if( map.empty() )  // a function without Variable: mark it as such
    map.push_back( Inf< Index >() );
  }

 // the ColVariable must be in the Block (or in its ancestors, as the
 // Variable of a LagBFunction are), and only OneVarConstraint are allowed;
 // checking the former would require scanning the whole tree, hence only
 // the type of the Variable of the Block is checked - - - - - - - - - - - -

 for( auto groups : { & f_Block->get_static_variable_groups() ,
		      & f_Block->get_dynamic_variable_groups() } )
  for( const auto & group : *groups )
   if( group && group->get_num_elements() &&
       ( ! group->elements_are< ColVariable >() ) )
    throw( std::logic_error( "SubGradSolver::set_Block: "
			     "some Variable is not a ColVariable" ) );

 for( auto groups : { & f_Block->get_static_constraint_groups() ,
		      & f_Block->get_dynamic_constraint_groups() } )
  for( const auto & group : *groups )
   if( group && group->get_num_elements() &&
       ( ! group->elements_are< OneVarConstraint >() ) )
    throw( std::logic_error( std::string( "SubGradSolver::set_Block: "
					  "unsupported type of " ) +
			     ( group->is_dynamic() ? "dynamic" : "static" ) +
			     " Constraint" ) );

 // the global pools: at least kPoolSz linearizations - - - - - - - - - - - -

 v_pool.assign( NrFi , false );
 v_dir_ok.assign( NrFi , false );
 v_slot.assign( NrFi , 0 );
 for( Index k = 0 ; k < NrFi ; ++k ) {
  auto f = v_c05f[ k ];
  const Index gps = f->get_int_par( C05Function::intGPMaxSz );
  const auto it = f_slot_base.find( f );
  if( ( it != f_slot_base.end() ) && ( it->second + kPoolSz <= gps ) ) {
   v_slot[ k ] = it->second;  // the names given before are still there
   v_pool[ k ] = true;
   continue;
   }
  try {
   f->set_par( C05Function::intGPMaxSz , int( gps + kPoolSz ) );
   v_slot[ k ] = gps;
   v_pool[ k ] = true;
   f_slot_base[ f ] = gps;
   }
  catch( std::exception & ) {}  // a C05Function without a global pool
  }

 // the bounds and the vectors - - - - - - - - - - - - - - - - - - - - - - -

 read_bounds();

 LambdaBar.resize( NumVar );
 for( Index i = 0 ; i < NumVar ; ++i )
  LambdaBar[ i ] = LamVcblr[ i ]->get_value();
 proj_fsb( LambdaBar );

 Lambda = LambdaBar;
 LambdaBest = LambdaBar;
 LambdaHat = LambdaBar;
 Gi.assign( NumVar , 0 );
 GiRaw.assign( NumVar , 0 );
 dir.assign( NumVar , 0 );
 dirRaw.assign( NumVar , 0 );

 Seq.clear();
 for( Index h = 0 ; h < NrFi + ( f_zero ? 1 : 0 ) ; ++h )
  Seq.push_back( h );

 }  // end( SubGradSolver::build_structure )

/*--------------------------------------------------------------------------*/

void SubGradSolver::clear_structure( void )
{
 v_c05f.clear();
 f_zero = nullptr;
 NrFi = NumVar = 0;
 LamVcblr.clear();
 v_loc2glob.clear();
 lb.clear();
 ub.clear();
 v_pool.clear();
 v_dir_ok.clear();
 v_slot.clear();
 LambdaBar.clear();
 Lambda.clear();
 LambdaBest.clear();
 LambdaHat.clear();
 Gi.clear();
 GiRaw.clear();
 dir.clear();
 dirRaw.clear();
 Seq.clear();
 f_eval_ok = false;
 f_need_rebuild = false;
 }

/*--------------------------------------------------------------------------*/

void SubGradSolver::read_bounds( void )
{
 lb.assign( NumVar , - INF );
 ub.assign( NumVar , INF );
 f_empty_box = false;

 for( Index i = 0 ; i < NumVar ; ++i ) {
  const auto var = LamVcblr[ i ];
  if( var->is_fixed() ) {
   lb[ i ] = ub[ i ] = var->get_value();
   continue;
   }

  lb[ i ] = var->get_lb();
  ub[ i ] = var->get_ub();

  for( Index j = 0 ; j < var->get_num_active() ; ++j ) {
   const auto c = dynamic_cast< const OneVarConstraint * >(
						      var->get_active( j ) );
   if( ( ! c ) || c->is_relaxed() )
    continue;
   lb[ i ] = std::max( lb[ i ] , VarValue( c->get_lhs() ) );
   ub[ i ] = std::min( ub[ i ] , VarValue( c->get_rhs() ) );
   }

  if( lb[ i ] > ub[ i ] )
   f_empty_box = true;
  }
 }  // end( SubGradSolver::read_bounds )

/*--------------------------------------------------------------------------*/

void SubGradSolver::build_rules( void )
{
 delete_rules();

 if( StepsizeName.empty() )
  throw( std::invalid_argument( "SubGradSolver::build_rules: "
				"strStepsize cannot be empty" ) );

 stepsize = Stepsize::new_Stepsize( StepsizeName , this );

 if( auto d = stepsize->deflection_role() ) {
  // the Stepsize is also the Deflection
  if( ( ! is_none( DeflectionName ) ) &&
      ( SMSpp_classname_normalise( std::string( DeflectionName ) ) !=
	stepsize->classname() ) ) {
   delete_rules();
   throw( std::invalid_argument( "SubGradSolver::build_rules: " +
				 StepsizeName + " is also the Deflection, "
				 "strDeflection cannot be " +
				 DeflectionName ) );
   }
  deflection = d;
  f_shared_rule = true;
  }
 else
  if( ! is_none( DeflectionName ) ) {
   deflection = Deflection::new_Deflection( DeflectionName , this );
   if( deflection->is_also_stepsize() ) {
    delete_rules();
    throw( std::invalid_argument( "SubGradSolver::build_rules: " +
				  DeflectionName + " is also a Stepsize, "
				  "it must be given in strStepsize" ) );
    }
   }

 f_rules_dirty = false;

 }  // end( SubGradSolver::build_rules )

/*--------------------------------------------------------------------------*/

void SubGradSolver::delete_rules( void )
{
 if( f_shared_rule )
  delete stepsize;
 else {
  delete stepsize;
  delete deflection;
  }

 stepsize = nullptr;
 deflection = nullptr;
 f_shared_rule = false;
 }

/*--------------------------------------------------------------------------*/

void SubGradSolver::reset_function_info( void )
{
 FiBar = FiLambda = FiBest = FiHat = INF;
 CertLB = - INF;
 Sigma = SigmaHat = 0;
 NrmGi = 0;
 MaxNrmGi = 0;
 DoSS = true;
 CSSCntr = CNSCntr = CSmallStep = 0;
 TotIter = 0;
 step = 0;
 std::fill( Gi.begin() , Gi.end() , 0 );
 std::fill( GiRaw.begin() , GiRaw.end() , 0 );
 GiCnst = 0;

 init_deflection();

 if( stepsize )   // if shared, this also formats the Deflection
  stepsize->format();

 }  // end( SubGradSolver::reset_function_info )

/*--------------------------------------------------------------------------*/

void SubGradSolver::init_deflection( void )
{
 NrmDir = dGk = dM1Gk = 0;
 NrmPrj = INF;
 Epsilon = HatEpsilon = 0;
 alpha = 1;
 dM1GkDone = false;
 f_dir_valid = false;
 ErgSum = 0;
 std::fill( dir.begin() , dir.end() , 0 );
 std::fill( dirRaw.begin() , dirRaw.end() , 0 );
 dirCnst = 0;
 v_dir_ok.assign( NrFi , false );

 if( deflection && ( ! f_shared_rule ) )
  deflection->format();

 }  // end( SubGradSolver::init_deflection )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR PARAMETERS --------------------------*/
/*--------------------------------------------------------------------------*/

void SubGradSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ):
   if( value < 0 )
    throw( std::invalid_argument(
			    "SubGradSolver::set_par: intMaxIter must be >= 0"
			    ) );
   MaxIter = value;
   return;
  case( intLogVerb ):  LogVerb = value; return;
  case( intEverykIt ): EverykIt = value; return;
  case( intProjection ):
   if( ( value < 0 ) || ( value > 7 ) )
    throw( std::invalid_argument(
		  "SubGradSolver::set_par: intProjection must be in 0 .. 7" )
	   );
   SGPar1 = value;
   return;
  case( intScheme ):
   if( ( value < 0 ) || ( value > 3 ) )
    throw( std::invalid_argument(
		      "SubGradSolver::set_par: intScheme must be in 0 .. 3" )
	   );
   SGPar3 = value;
   return;
  case( intKeepHat ): SGPar4 = value; return;
  case( intSeed ):
   SGPar5 = value;
   myrandom.seed( SGPar5 );
   return;
  case( intResetAlg ): RstAlg = value; return;
  case( intPrimalRec ):
   if( ( value < 0 ) || ( value > 1 ) )
    throw( std::invalid_argument(
		    "SubGradSolver::set_par: intPrimalRec must be 0 or 1" ) );
   if( value != PrimalRec ) {  // the combination kept so far is another
    PrimalRec = value;
    ErgSum = 0;
    v_dir_ok.assign( NrFi , false );
    }
   return;
  default:
   if( ( par >= intColorGreen ) && ( par < intLastParSGS ) ) {
    if( value < 0 )
     throw( std::invalid_argument( "SubGradSolver::set_par: " +
				   int_par_idx2str( par ) + " must be >= 0" )
	    );
    v_rule_int[ par - intColorGreen ] = value;
    return;
    }
   CDASolver::set_par( par , value );
  }
 }  // end( SubGradSolver::set_par( int ) )

/*--------------------------------------------------------------------------*/

void SubGradSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblMaxTime ):  MaxTime = value; return;
  case( dblEveryTTm ): EveryTTm = value; return;
  case( dblRelAcc ):
   if( value < 0 )
    throw( std::invalid_argument(
			     "SubGradSolver::set_par: dblRelAcc must be >= 0"
			     ) );
   RelAcc = value;
   return;
  case( dblAbsAcc ):
   if( value < 0 )
    throw( std::invalid_argument(
			     "SubGradSolver::set_par: dblAbsAcc must be >= 0"
			     ) );
   AbsAcc = value;
   return;
  case( dblUpCutOff ): UpCutOff = value; return;
  case( dblLwCutOff ): LwCutOff = value; return;
  case( dbltStar ):
   if( value <= 0 )
    throw( std::invalid_argument(
			       "SubGradSolver::set_par: dbltStar must be > 0"
			       ) );
   tStar = value;
   return;
  case( dblNZEps ):
   if( value < 0 )
    throw( std::invalid_argument(
			       "SubGradSolver::set_par: dblNZEps must be >= 0"
			       ) );
   NZEps = value;
   return;
  case( dblIncrFactor ):
   SGPar2 = std::max( value , double( 0 ) );
   return;
  default:
   if( ( par >= dblLipFactor ) && ( par < dblLastParSGS ) ) {
    if( value < 0 )
     throw( std::invalid_argument( "SubGradSolver::set_par: " +
				   dbl_par_idx2str( par ) + " must be >= 0" )
	    );
    v_rule_dbl[ par - dblLipFactor ] = value;
    return;
    }
   CDASolver::set_par( par , value );
  }
 }  // end( SubGradSolver::set_par( double ) )

/*--------------------------------------------------------------------------*/

void SubGradSolver::set_par( idx_type par , std::string && value )
{
 switch( par ) {
  case( strStepsize ):
   if( value != StepsizeName ) {
    StepsizeName = std::move( value );
    f_rules_dirty = true;
    }
   return;
  case( strDeflection ):
   if( value != DeflectionName ) {
    DeflectionName = std::move( value );
    f_rules_dirty = true;
    }
   return;
  default:
   CDASolver::set_par( par , std::move( value ) );
  }
 }  // end( SubGradSolver::set_par( std::string ) )

/*--------------------------------------------------------------------------*/

int SubGradSolver::get_dflt_int_par( idx_type par ) const
{
 if( ( par >= intProjection ) && ( par < intLastParSGS ) )
  return( SGS_dflt_int_par[ par - intProjection ] );

 return( CDASolver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double SubGradSolver::get_dflt_dbl_par( idx_type par ) const
{
 if( ( par >= dbltStar ) && ( par < dblLastParSGS ) )
  return( SGS_dflt_dbl_par[ par - dbltStar ] );

 return( CDASolver::get_dflt_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & SubGradSolver::get_dflt_str_par( idx_type par ) const
{
 if( ( par >= strStepsize ) && ( par < strLastParSGS ) )
  return( SGS_dflt_str_par[ par - strStepsize ] );

 return( CDASolver::get_dflt_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int SubGradSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ):    return( MaxIter );
  case( intLogVerb ):    return( LogVerb );
  case( intEverykIt ):   return( EverykIt );
  case( intProjection ): return( SGPar1 );
  case( intScheme ):     return( SGPar3 );
  case( intKeepHat ):    return( SGPar4 );
  case( intSeed ):       return( SGPar5 );
  case( intResetAlg ):   return( RstAlg );
  case( intPrimalRec ):  return( PrimalRec );
  default:
   if( ( par >= intColorGreen ) && ( par < intLastParSGS ) )
    return( v_rule_int[ par - intColorGreen ] );
  }

 return( CDASolver::get_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

double SubGradSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ):    return( MaxTime );
  case( dblEveryTTm ):   return( EveryTTm );
  case( dblRelAcc ):     return( RelAcc );
  case( dblAbsAcc ):     return( AbsAcc );
  case( dblUpCutOff ):   return( UpCutOff );
  case( dblLwCutOff ):   return( LwCutOff );
  case( dbltStar ):      return( tStar );
  case( dblNZEps ):      return( NZEps );
  case( dblIncrFactor ): return( SGPar2 );
  default:
   if( ( par >= dblLipFactor ) && ( par < dblLastParSGS ) )
    return( v_rule_dbl[ par - dblLipFactor ] );
  }

 return( CDASolver::get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & SubGradSolver::get_str_par( idx_type par ) const
{
 if( par == strStepsize )
  return( StepsizeName );

 if( par == strDeflection )
  return( DeflectionName );

 return( CDASolver::get_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

SubGradSolver::idx_type SubGradSolver::int_par_str2idx(
					     const std::string & name ) const
{
 const auto it = SGS_int_pars_map.find( name );
 if( it != SGS_int_pars_map.end() )
  return( it->second );

 return( CDASolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

SubGradSolver::idx_type SubGradSolver::dbl_par_str2idx(
					     const std::string & name ) const
{
 const auto it = SGS_dbl_pars_map.find( name );
 if( it != SGS_dbl_pars_map.end() )
  return( it->second );

 return( CDASolver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

SubGradSolver::idx_type SubGradSolver::str_par_str2idx(
					     const std::string & name ) const
{
 const auto it = SGS_str_pars_map.find( name );
 if( it != SGS_str_pars_map.end() )
  return( it->second );

 return( CDASolver::str_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & SubGradSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intProjection ) && ( idx < intLastParSGS ) )
  return( SGS_int_pars_str[ idx - intProjection ] );

 return( CDASolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & SubGradSolver::dbl_par_idx2str( idx_type idx ) const
{
 if( ( idx >= dbltStar ) && ( idx < dblLastParSGS ) )
  return( SGS_dbl_pars_str[ idx - dbltStar ] );

 return( CDASolver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & SubGradSolver::str_par_idx2str( idx_type idx ) const
{
 if( ( idx >= strStepsize ) && ( idx < strLastParSGS ) )
  return( SGS_str_pars_str[ idx - strStepsize ] );

 return( CDASolver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

int SubGradSolver::compute( bool changedvars )
{
 lock();  // ensure no concurrent accesses

 if( Result == kStillRunning ) {
  unlock();
  throw( std::logic_error( "SubGradSolver::compute: called within itself" )
	 );
  }

 if( ! f_Block ) {
  Result = kBlockLocked;
  unlock();
  return( Result );
  }

 c_start = std::chrono::system_clock::now();

 // process the outstanding Modification - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_mod.empty() ) {
  bool owned = f_Block->is_owned_by( f_id );
  if( ( ! owned ) && ( ! f_Block->read_lock() ) ) {
   Result = kBlockLocked;
   unlock();
   return( Result );
   }

  try {
   process_outstanding_Modification();
   }
  catch( ... ) {
   if( ! owned )
    f_Block->read_unlock();
   unlock();
   throw;
   }

  if( ! owned )
   f_Block->read_unlock();
  }

 // the rules, if changed - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_rules_dirty || ( ! stepsize ) ) {
  try {
   build_rules();
   }
  catch( ... ) {
   unlock();
   throw;
   }
  // the function has not changed, what the old rules knew of it has
  init_deflection();
  stepsize->format();
  TotIter = 0;
  }

 ++SCalls;
 ParIter = 0;

 if( MaxIter == 0 ) {
  Result = kStopIter;
  unlock();
  return( Result );
  }

 if( f_empty_box ) {
  SGLOG( 0 , std::endl << "SubGradSolver: empty feasible set" << std::endl
	 );
  Result = kInfeasible;
  unlock();
  return( Result );
  }

 // reset what has to be at each call- - - - - - - - - - - - - - - - - - - -

 if( SCalls > 1 ) {
  if( RstAlg & RstCenter ) {
   for( Index i = 0 ; i < NumVar ; ++i )
    LambdaBar[ i ] = LamVcblr[ i ]->get_value();
   proj_fsb( LambdaBar );
   FiBar = INF;
   }

  if( RstAlg & RstRules )
   reset_function_info();
  }

 // the scheme actually used and the incremental variant - - - - - - - - - -

 f_scheme = SGPar3;
 if( ! deflection )
  f_scheme |= StepRST;
 if( f_shared_rule )  // the deflection makes the step, which has no beta
  f_scheme = StepRST;

 const Index NFun = NrFi + ( f_zero ? 1 : 0 );
 NItIncr = Index( std::ceil( double( NFun ) * SGPar2 ) );
 if( NItIncr && deflection ) {
  unlock();
  throw( std::logic_error( "SubGradSolver::compute: the incremental "
			   "variant cannot be deflected" ) );
  }
 InnIter = false;

 // the lower bounds of the components, which do not change during the
 // call, and those of the Block
 VarValue CompLB = - INF;
 if( ! f_zero ) {
  CompLB = 0;
  for( auto f : v_c05f ) {
   CompLB += f_convex ? f->get_global_lower_bound()
                      : - f->get_global_upper_bound();
   if( CompLB == - INF )
    break;
   }
  }

 auto updt_bounds = [ & ]( void ) {
  TrueLB = f_convex ? f_Block->get_valid_lower_bound( false )
                    : - f_Block->get_valid_upper_bound( false );
  TrueLB = std::max( TrueLB , CompLB );
  CondLB = f_convex ? f_Block->get_valid_lower_bound( true )
                    : - f_Block->get_valid_upper_bound( true );
  };

 updt_bounds();

 // the cutoffs, for the minimized function: a point at least as good as
 // PntCutOff, or a lower bound at least as large as BndCutOff, stop
 const VarValue PntCutOff = f_convex ? LwCutOff : - UpCutOff;
 const VarValue BndCutOff = f_convex ? UpCutOff : - LwCutOff;

 Result = kStillRunning;
 double lastETT = 0;

 SGLOG( 0 , std::endl << "SubGradSolver: " << NumVar << " variables, "
	<< NrFi << ( f_zero ? " + 1" : "" ) << " components, "
	<< ( stepsize ? stepsize->classname() : "?" ) << " + "
	<< ( deflection ? deflection->classname() : "no deflection" )
	<< ( NItIncr ? ", " + std::to_string( NItIncr ) + " incremental" :
	     "" ) << std::endl );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( ; ; ) {
  // check the time - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( MaxTime < INF ) && ( get_elapsed_time() > MaxTime ) ) {
   Result = kStopTime;
   break;
   }

  // time-periodic events - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( EveryTTm && ( get_elapsed_time() >= lastETT + EveryTTm ) ) {
   for( auto & ev : v_events[ eEveryTTime ] ) {
    auto res = ev();
    if( res == eStopOK ) { Result = kOK; break; }
    if( res == eStopError ) { Result = kError; break; }
    }
   lastETT = get_elapsed_time();
   if( Result != kStillRunning )
    break;
   }

  // the direction and the stepsize, if the center is known - - - - - - - - -

  if( FiBar < INF ) {
   form_D();
   save_dir();

   // the rule needs the scalar product of the previous direction with the
   // next subgradient, and without Deflection the direction is Gi
   if( ( ! deflection ) && stepsize->needs_dM1Gk() )
    dir = Gi;

   // a stop condition regarding consecutive "short" stepsizes
   if( step <= SmallStep * tStar ) {
    if( ++CSmallStep > ( NItIncr ? MaxSmallStep * NFun : MaxSmallStep ) ) {
     SGLOG( 0 , "SubGradSolver: the stepsize stays tiny" << std::endl );
     Result = kLowPrecision;
     break;
     }
    }
   else
    CSmallStep = 0;
   }

  // check for optimality - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ! InnIter ) {
   const bool cert = ( FiBest < INF ) &&
                     ( FiBest - global_lb() <= max_error( FiBest ) );
   if( cert || ( ( FiBar < INF ) && is_optimal() ) ) {
    int res = eContinue;
    for( auto & ev : v_events[ eBeforeTermination ] )
     if( ( res = ev() ) != eContinue )
      break;

    if( res == eStopError ) {
     Result = kError;
     break;
     }

    if( res != eForceContinue ) {
     SGLOG( 0 , "SubGradSolver: optimal"
	    << ( cert ? " (certified)" : "" ) << std::endl );
     Result = kOK;
     break;
     }
    }
   }

  // iteration-periodic events- - - - - - - - - - - - - - - - - - - - - - - -

  if( EverykIt && ParIter && ( ! ( ParIter % EverykIt ) ) ) {
   for( auto & ev : v_events[ eEverykIteration ] ) {
    auto res = ev();
    if( res == eStopOK ) { Result = kOK; break; }
    if( res == eStopError ) { Result = kError; break; }
    }
   if( Result != kStillRunning )
    break;
   }

  // the iterations limit - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( Index( MaxIter ) <= ParIter ) {
   Result = kStopIter;
   break;
   }

  // the tentative point (not projected yet)- - - - - - - - - - - - - - - - -

  if( FiBar < INF )
   form_lambda1();
  else
   Lambda = LambdaBar;

  ++ParIter;
  ++TotIter;

  // is it an inner iteration, and of which function? - - - - - - - - - - - -

  InnIter = NItIncr ? ( GiEvaltns % ( NItIncr + 1 ) ) != 0 : false;
  Index comp = Inf< Index >();
  if( InnIter ) {
   comp = ( ( GiEvaltns % ( NItIncr + 1 ) ) - 1 ) % NFun;
   if( ! comp )
    std::shuffle( Seq.begin() , Seq.end() , myrandom );
   comp = Seq[ comp ];
   }

  // project the point on the box and evaluate the function there - - - - - -

  LHasProj = proj_fsb( Lambda );

  if( ! Fi_and_Gi( comp ) ) {
   Result = kError;
   break;
   }

  if( FiLambda == - INF ) {
   SGLOG( 0 , "SubGradSolver: the function is unbounded below" << std::endl
	  );
   Result = kUnbounded;
   break;
   }

  // project the subgradient over the tangent cone at the point - - - - - - -

  if( SGPar1 & PrjSG )
   proj_tng( Gi , Lambda );

  NrmGi = norm2( Gi );

  if( stepsize->needs_dM1Gk() ) {
   dM1Gk = dot( dir , Gi );
   dM1GkDone = true;
   }
  else
   dM1GkDone = false;

  // the iteration log- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( f_log && ( LogVerb > 1 ) ) {
   *f_log << std::setw( 6 ) << ParIter << ( InnIter ? "i" : " " );
   sci( *f_log , 8 );
   *f_log << " Fi " << FiLambda << " FiBest " << FiBest;
   sci( *f_log , 2 );
   *f_log << " step " << step << " alpha " << alpha
	  << " |g| " << std::sqrt( NrmGi ) << " |d| " << std::sqrt( NrmDir )
	  << " eps " << Epsilon;
   if( global_lb() > - INF ) {
    sci( *f_log , 8 );
    *f_log << " LB " << global_lb();
    }
   *f_log << std::endl;
   }

  // the bounds, and the certified one when the subgradient is complete - - -

  if( ! InnIter ) {
   updt_bounds();
   if( ! deflection )
    update_cert_bound();
   }

  // move the stability center or not - - - - - - - - - - - - - - - - - - - -

  if( ( FiBar == INF ) || ( ! deflection ) ) {
   goto_lambda1();
   DoSS = true;
   }
  else
   if( ( DoSS = deflection->do_SS() ) ) {  // serious step
    goto_lambda1();
    ++CSSCntr;
    CNSCntr = 0;
    SGLOG( 2 , "        SS[" << CSSCntr << "]" << std::endl );
    }
   else {                                   // null step
    update_sigma();
    ++CNSCntr;
    CSSCntr = 0;
    SGLOG( 2 , "        NS[" << CNSCntr << "]" << std::endl );
    }

  // the conditional lower bound and the cutoffs- - - - - - - - - - - - - - -

  if( ( CondLB > - INF ) && ( FiBest < CondLB - max_error( CondLB ) ) ) {
   SGLOG( 0 , "SubGradSolver: below the conditional lower bound"
	  << std::endl );
   Result = kUnbounded;
   break;
   }

  if( ( FiBest <= PntCutOff ) || ( global_lb() >= BndCutOff ) ) {
   SGLOG( 0 , "SubGradSolver: cutoff" << std::endl );
   Result = kCutOff;
   break;
   }
  }  // end( main cycle )

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_log && ( LogVerb > 0 ) ) {
  *f_log << "SubGradSolver: status " << Result << " after " << ParIter
	 << " iterations (" << FiEvaltns << " evaluations), ";
  sci( *f_log , 10 );
  *f_log << "best " << ( f_convex ? FiBest : - FiBest ) << ", bound "
	 << ( f_convex ? global_lb() : - global_lb() );
  *f_log << std::setprecision( 4 ) << std::fixed << " [" << get_elapsed_time()
	 << " s]" << std::endl;
  }

 unlock();
 return( Result );

 }  // end( SubGradSolver::compute )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue SubGradSolver::get_lb( void )
{
 if( f_convex )
  return( global_lb() );

 const auto v = sol_value();
 return( v < INF ? - v : - INF );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue SubGradSolver::get_ub( void )
{
 if( f_convex )
  return( sol_value() );

 return( - global_lb() );
 }

/*--------------------------------------------------------------------------*/

bool SubGradSolver::has_var_solution( void )
{
 return( f_Block != nullptr );
 }

/*--------------------------------------------------------------------------*/

bool SubGradSolver::is_var_feasible( void )
{
 return( ( sol_value() < INF ) && ( ! f_empty_box ) );
 }

/*--------------------------------------------------------------------------*/

void SubGradSolver::get_var_solution( Configuration * solc )
{
 // a pending change of the ColVariable may have deleted some of them
 if( ! v_mod.empty() )
  process_outstanding_Modification();

 lambda_to_vars( FiBest < INF ? LambdaBest : LambdaBar );

 }  // end( SubGradSolver::get_var_solution )

/*--------------------------------------------------------------------------*/

bool SubGradSolver::has_dual_solution( void )
{
 if( ! f_eval_ok )
  return( false );

 for( Index k = 0 ; k < NrFi ; ++k ) {
  if( ! v_pool[ k ] )
   return( false );
  const Index slot = v_slot[ k ] +
                     ( v_dir_ok[ k ] ? DirPos : kSlotG );
  if( ! v_c05f[ k ]->is_linearization_there( slot ) )
   return( false );
  }

 return( true );
 }

/*--------------------------------------------------------------------------*/

bool SubGradSolver::is_dual_feasible( void )
{
 return( global_lb() > - INF );
 }

/*--------------------------------------------------------------------------*/

void SubGradSolver::get_dual_solution( Configuration * solc )
{
 auto set_k = [ this ]( Index k ) {
  if( ! v_pool[ k ] )
   throw( std::logic_error( "SubGradSolver::get_dual_solution: component "
			    + std::to_string( k ) + " has no global pool" ) );
  const Index slot = v_slot[ k ] +
                     ( v_dir_ok[ k ] ? DirPos : kSlotG );
  v_c05f[ k ]->set_important_linearization(
			   LinearCombination( { std::pair( slot , 1.0 ) } ) );
  };

 if( auto c = dynamic_cast< SimpleConfiguration< std::vector< int > > * >(
								    solc ) ) {
  for( auto k : c->f_value ) {
   if( ( k < 0 ) || ( Index( k ) >= NrFi ) )
    throw( std::invalid_argument( "SubGradSolver::get_dual_solution: "
				  "invalid index " + std::to_string( k ) ) );
   set_k( k );
   }
  return;
  }

 for( Index k = 0 ; k < NrFi ; ++k )
  set_k( k );

 }  // end( SubGradSolver::get_dual_solution )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

bool SubGradSolver::process_outstanding_Modification( void )
{
 // atomically move v_mod in a flattened temporary list - - - - - - - - - - -

 Lst_sp_Mod v_mod_tmp;

 std::function< void( const sp_Mod & ) > flatten =
  [ & ]( const sp_Mod & mod ) {
   if( auto g = std::dynamic_pointer_cast< GroupModification >( mod ) )
    for( const auto & sm : g->sub_Modifications() )
     flatten( sm );
   else
    v_mod_tmp.push_back( mod );
   };

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;  // try to acquire lock, spin on failure

 for( const auto & mod : v_mod )
  flatten( mod );

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

 // classify them: a change of the set of ColVariable or of their bounds
 // requires to read the structure again, a change of a function makes all
 // that is known about the function useless- - - - - - - - - - - - - - - -

 bool rebuild = false;
 bool funchg = false;

 auto ours = [ this ]( const Function * f ) -> Index {
  for( Index k = 0 ; k < NrFi ; ++k )
   if( v_c05f[ k ] == f )
    return( k );
  return( ( f && ( f == f_zero ) ) ? NrFi : Inf< Index >() );
  };

 for( const auto & mod : v_mod_tmp ) {
  if( std::dynamic_pointer_cast< NBModification >( mod ) ) {
   rebuild = true;
   continue;
   }

  if( std::dynamic_pointer_cast< FunctionModVars >( mod ) ) {
   rebuild = true;
   continue;
   }

  if( auto cmod = std::dynamic_pointer_cast< C05FunctionMod >( mod ) ) {
   const auto k = ours( cmod->function() );
   if( k == Inf< Index >() )
    continue;
   switch( cmod->type() ) {
    case( C05FunctionMod::NothingChanged ):
    case( C05FunctionMod::GlobalPoolAdded ):
     break;
    case( C05FunctionMod::GlobalPoolRemoved ):
     // the pool may have lost the direction
     if( k < NrFi )
      v_dir_ok[ k ] = false;
     break;
    default:
     funchg = true;
    }
   if( ! std::isnan( cmod->shift() ) && ( cmod->shift() != 0 ) )
    funchg = true;
   continue;
   }

  if( auto fmod = std::dynamic_pointer_cast< FunctionMod >( mod ) ) {
   if( ours( fmod->function() ) != Inf< Index >() )
    funchg = true;
   continue;
   }

  if( std::dynamic_pointer_cast< BlockModAD >( mod ) ||
      std::dynamic_pointer_cast< BlockMod >( mod ) ||
      std::dynamic_pointer_cast< VariableMod >( mod ) ||
      std::dynamic_pointer_cast< ConstraintMod >( mod ) ) {
   rebuild = true;
   continue;
   }

  // anything else is of no concern
  }

 if( rebuild ) {
  // keep the center by ColVariable: the old pointers are only used as keys,
  // as some of the ColVariable may not exist any longer
  std::unordered_map< const ColVariable * , VarValue > old;
  for( Index i = 0 ; i < NumVar ; ++i )
   old.emplace( LamVcblr[ i ] , LambdaBar[ i ] );

  clear_structure();
  build_structure( false );

  for( Index i = 0 ; i < NumVar ; ++i ) {
   auto it = old.find( LamVcblr[ i ] );
   if( it != old.end() )
    LambdaBar[ i ] = it->second;
   }
  proj_fsb( LambdaBar );
  Lambda = LambdaBest = LambdaHat = LambdaBar;

  reset_function_info();
  return( true );
  }

 if( funchg )
  reset_function_info();

 return( funchg );

 }  // end( SubGradSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/

void SubGradSolver::lambda_to_vars( const Vec_VarValue & L )
{
 for( Index i = 0 ; i < NumVar ; ++i )
  LamVcblr[ i ]->set_value( L[ i ] );
 }

/*--------------------------------------------------------------------------*/

bool SubGradSolver::eval_function( Index h , VarValue & val ,
				   VarValue & cnst )
{
 const bool zero = ( h >= NrFi );
 auto f = zero ? f_zero : v_c05f[ h ];

 if( ( ! zero ) && ( MaxTime < INF ) )
  f->set_par( dblMaxTime , std::max( MaxTime - get_elapsed_time() , 0.0 ) );

 const int st = f->compute( true );
 if( ( st <= kUnEval ) || ( ( st >= kError ) && ( st != kLowPrecision ) ) ) {
  SGLOG( 0 , "SubGradSolver: error " << st << " evaluating component "
	 << h << std::endl );
  return( false );
  }

 val = f_convex ? f->get_upper_estimate() : - f->get_lower_estimate();

 if( val == INF ) {
  // the point is out of the domain of the function: only a vertical
  // linearization could be given, which SubGradSolver cannot handle
  SGLOG( 0 , "SubGradSolver: component " << h << " is + INF, vertical "
	 "linearizations are not handled" << std::endl );
  return( false );
  }

 if( val == - INF )  // the function is unbounded below
  return( true );

 if( ! f->has_linearization( true ) ) {
  SGLOG( 0 , "SubGradSolver: component " << h << " has no diagonal "
	 "linearization" << std::endl );
  return( false );
  }

 const auto & map = v_loc2glob[ zero ? NrFi : h ];
 const Index nloc = f->get_num_active_var();
 v_buf.resize( nloc );
 if( nloc )
  f->get_linearization_coefficients( v_buf.data() );

 if( map.empty() ) {
  if( f_convex )
   for( Index i = 0 ; i < nloc ; ++i )
    GiRaw[ i ] += v_buf[ i ];
  else
   for( Index i = 0 ; i < nloc ; ++i )
    GiRaw[ i ] -= v_buf[ i ];
  }
 else
  if( nloc ) {
   if( f_convex )
    for( Index i = 0 ; i < nloc ; ++i )
     GiRaw[ map[ i ] ] += v_buf[ i ];
   else
    for( Index i = 0 ; i < nloc ; ++i )
     GiRaw[ map[ i ] ] -= v_buf[ i ];
   }

 cnst = rs( f->get_linearization_constant() );

 if( ( ! zero ) && v_pool[ h ] )
  f->store_linearization( v_slot[ h ] + kSlotG );

 return( true );

 }  // end( SubGradSolver::eval_function )

/*--------------------------------------------------------------------------*/

bool SubGradSolver::Fi_and_Gi( Index wFi )
{
 lambda_to_vars( Lambda );

 std::fill( GiRaw.begin() , GiRaw.end() , 0 );
 GiCnst = 0;

 if( wFi == Inf< Index >() ) {  // all the functions- - - - - - - - - - - - -
  VarValue tot = 0;
  const Index NFun = NrFi + ( f_zero ? 1 : 0 );
  for( Index h = 0 ; h < NFun ; ++h ) {
   VarValue val = 0;
   VarValue cnst = 0;
   if( ! eval_function( h , val , cnst ) )
    return( false );
   if( val == - INF ) {
    FiLambda = - INF;
    return( true );
    }
   tot += val;
   GiCnst += cnst;
   }

  FiLambda = tot;
  ++FiEvaltns;
  f_eval_ok = true;

  if( FiLambda < FiBest ) {
   FiBest = FiLambda;
   LambdaBest = Lambda;
   }

  Sigma = FiLambda - ( GiCnst + dot( GiRaw , Lambda ) );
  MaxNrmGi = std::max( MaxNrmGi , norm2( GiRaw ) );
  }
 else {                       // only the function wFi- - - - - - - - - - - -
  VarValue val = 0;
  if( ! eval_function( wFi , val , GiCnst ) )
   return( false );
  if( val == - INF ) {
   FiLambda = - INF;
   return( true );
   }
  Sigma = val - ( GiCnst + dot( GiRaw , Lambda ) );
  }

 if( SGPar4 )
  SigmaHat = Sigma;

 Gi = GiRaw;
 ++GiEvaltns;

 return( true );

 }  // end( SubGradSolver::Fi_and_Gi )

/*--------------------------------------------------------------------------*/

void SubGradSolver::form_D( void )
{
 // after a null step the subgradient has been projected at the tentative
 // point, while it must be at the stability center
 if( ( ! DoSS ) && ( SGPar1 & PrjSG ) ) {
  Gi = GiRaw;
  proj_tng( Gi , LambdaBar );
  NrmGi = norm2( Gi );
  dM1GkDone = false;
  }

 // project the previous direction, if so required
 if( deflection && ( SGPar1 & PrjDIRM1 ) ) {
  if( SGPar1 & PrjDIR ) {  // projected already, at the previous center
   if( DoSS ) {            // which has changed
    proj_tng( dir , LambdaBar );
    NrmDir = norm2( dir );
    dM1GkDone = false;
    }
   }
  else {
   proj_tng( dir , LambdaBar );
   NrmDir = norm2( dir );
   dM1GkDone = false;
   }
  }

 // the scalar product between the previous direction and the subgradient
 if( deflection ) {
  if( ! dM1GkDone )
   dM1Gk = dot( Gi , dir );
  }
 else {  // without Deflection the direction is the subgradient
  dGk = NrmDir = NrmGi;
  if( ! dM1GkDone )
   dM1Gk = dGk;
  }
 dM1GkDone = true;

 // the stepsize and the deflection- - - - - - - - - - - - - - - - - - - - -

 if( f_scheme & StepRST ) {  // stepsize-restricted: deflection first
  if( deflection ) {
   deflection->new_deflection();
   alpha = deflection->get_deflection();
   }
  }
 else {                      // deflection-restricted: stepsize first
  dGk = dM1Gk;
  if( ! InnIter )
   stepsize->new_step();
  step = stepsize->get_stepsize( InnIter );

  deflection->new_deflection();
  alpha = deflection->get_deflection();

  // the safe rule
  if( ( f_scheme & SafeRULE ) && ( ! InnIter ) ) {
   const VarValue den = FiBar - stepsize->get_level() + step * NrmDir;
   if( den > 0 ) {
    const VarValue target = ( step * NrmDir ) / den;
    alpha = std::min( VarValue( 1 ) , std::max( alpha , target ) );
    }
   }
  }

 // with no valid previous direction, the direction is the subgradient
 if( ! f_dir_valid )
  alpha = 1;

 }  // end( SubGradSolver::form_D )

/*--------------------------------------------------------------------------*/

void SubGradSolver::save_dir( void )
{
 // the linearization error of the new direction at the center
 if( alpha == 1 )
  Epsilon = Sigma;
 else
  Epsilon = alpha * Sigma + ( 1 - alpha ) * Epsilon;

 if( SGPar4 && ( ! NItIncr ) )
  chg_lambda_hat();

 if( deflection ) {
  // the direction is the convex combination of the old one and Gi, and so
  // is the aggregate linearization
  if( alpha == 1 ) {
   dir = Gi;
   dirRaw = GiRaw;
   dirCnst = GiCnst;
   }
  else {
   comb( dir , 1 - alpha , Gi , alpha );
   comb( dirRaw , 1 - alpha , GiRaw , alpha );
   dirCnst = ( 1 - alpha ) * dirCnst + alpha * GiCnst;
   }

  f_dir_valid = true;
  update_cert_bound();

  if( SGPar1 & PrjDIR )
   proj_tng( dir , LambdaBar );

  NrmDir = norm2( dir );
  dGk = dot( Gi , dir );
  }

 // the norm of the direction projected on the tangent cone at the center,
 // whatever is projected in the iterations: the direction is an
 // Epsilon-subgradient of the function plus the indicator of the box, hence
 // if the projection is "0" the center is Epsilon-optimal
 const auto & d = deflection ? dir : Gi;
 NrmPrj = 0;
 for( Index i = 0 ; i < NumVar ; ++i )
  if( ! ( ( ( LambdaBar[ i ] >= ub[ i ] ) && ( d[ i ] <= 0 ) ) ||
	  ( ( LambdaBar[ i ] <= lb[ i ] ) && ( d[ i ] >= 0 ) ) ) )
   NrmPrj += d[ i ] * d[ i ];

 if( ( ! InnIter ) && ( std::sqrt( NrmPrj ) <= NZEps ) &&
     ( FiBar - Epsilon > CertLB ) )
  CertLB = FiBar - Epsilon;

 // the stepsize in the stepsize-restricted scheme- - - - - - - - - - - - - -

 if( f_scheme & StepRST ) {
  if( ! InnIter )
   stepsize->new_step();

  step = stepsize->get_stepsize( InnIter );

  if( f_scheme & SafeRULE ) {
   stepsize->set_max_beta( alpha );
   step = stepsize->get_stepsize( InnIter );
   }
  }

 // the convex combination of linearizations for the primal recovery, only
 // at outer iterations, whose subgradient is that of all the components
 if( ! InnIter )
  store_primal();

 }  // end( SubGradSolver::save_dir )

/*--------------------------------------------------------------------------*/

void SubGradSolver::store_primal( void )
{
 // the weight of the latest linearization: the deflection coefficient, or
 // that of the stepsize in the ergodic average
 VarValue w = 1;
 if( deflection && ( ! PrimalRec ) )
  w = alpha;
 else {
  if( step <= 0 )  // no move, no weight
   return;
  ErgSum += step;
  w = step / ErgSum;
  }

 const Index newpos = ( DirPos == 1 ) ? 2 : 1;

 for( Index k = 0 ; k < NrFi ; ++k ) {
  if( ! v_pool[ k ] )
   continue;

  auto f = v_c05f[ k ];
  const Index g = v_slot[ k ] + kSlotG;
  const Index d = v_slot[ k ] + DirPos;
  if( ! f->is_linearization_there( g ) ) {
   v_dir_ok[ k ] = false;
   continue;
   }

  if( w == 1 ) {
   f->store_combination_of_linearizations(
		      LinearCombination( { std::pair( g , 1.0 ) } ) ,
		      v_slot[ k ] + newpos );
   v_dir_ok[ k ] = true;
   }
  else
   if( v_dir_ok[ k ] && f->is_linearization_there( d ) )
    f->store_combination_of_linearizations(
	      LinearCombination( { std::pair( g , w ) ,
				   std::pair( d , 1 - w ) } ) ,
	      v_slot[ k ] + newpos );
   else  // the old combination is lost, and so is the new one
    v_dir_ok[ k ] = false;
  }

 DirPos = newpos;

 }  // end( SubGradSolver::store_primal )

/*--------------------------------------------------------------------------*/

void SubGradSolver::goto_lambda1( void )
{
 // the linearization error of the direction at the new center
 if( ( FiBar < INF ) && deflection ) {
  VarValue SclPrd = 0;
  if( ! LHasProj )
   SclPrd = step * NrmDir;
  else
   for( Index i = 0 ; i < NumVar ; ++i )
    SclPrd += dir[ i ] * ( LambdaBar[ i ] - Lambda[ i ] );

  Epsilon += FiLambda - FiBar + SclPrd;
  }

 LambdaBar = Lambda;
 if( ! InnIter )
  FiBar = FiLambda;

 }  // end( SubGradSolver::goto_lambda1 )

/*--------------------------------------------------------------------------*/

void SubGradSolver::form_lambda1( void )
{
 const auto & d = deflection ? dir : Gi;
 for( Index i = 0 ; i < NumVar ; ++i )
  Lambda[ i ] = LambdaBar[ i ] - step * d[ i ];
 }

/*--------------------------------------------------------------------------*/

bool SubGradSolver::proj_fsb( Vec_VarValue & Lpt )
{
 bool chgd = false;
 for( Index i = 0 ; i < NumVar ; ++i )
  if( Lpt[ i ] > ub[ i ] ) {
   Lpt[ i ] = ub[ i ];
   chgd = true;
   }
  else
   if( Lpt[ i ] < lb[ i ] ) {
    Lpt[ i ] = lb[ i ];
    chgd = true;
    }

 return( chgd );
 }

/*--------------------------------------------------------------------------*/

void SubGradSolver::proj_tng( Vec_VarValue & Gpt , const Vec_VarValue & Lpt )
{
 // the step is along - Gpt: a component that would leave the box from the
 // face where the point is is zeroed
 for( Index j = 0 ; j < NumVar ; ++j )
  if( ( ( Lpt[ j ] >= ub[ j ] ) && ( Gpt[ j ] <= 0 ) ) ||
      ( ( Lpt[ j ] <= lb[ j ] ) && ( Gpt[ j ] >= 0 ) ) )
   Gpt[ j ] = 0;
 }

/*--------------------------------------------------------------------------*/

void SubGradSolver::update_sigma( void )
{
 // the linearization error of Gi at the center
 Sigma += FiBar - FiLambda;
 if( ! LHasProj )
  Sigma -= step * dot( Gi , dir );
 else
  for( Index i = 0 ; i < NumVar ; ++i )
   Sigma += Gi[ i ] * ( Lambda[ i ] - LambdaBar[ i ] );
 }

/*--------------------------------------------------------------------------*/

void SubGradSolver::chg_lambda_hat( void )
{
 if( alpha == 1 ) {
  HatEpsilon = SigmaHat;
  LambdaHat = Lambda;
  FiHat = FiLambda;
  return;
  }

 // the linearization error of the previous direction at LambdaHat
 VarValue SclPrd = 0;
 for( Index i = 0 ; i < NumVar ; ++i )
  SclPrd += dir[ i ] * ( LambdaHat[ i ] - Lambda[ i ] );

 HatEpsilon += alpha * ( FiLambda - FiHat + SclPrd );

 // the linearization error of the subgradient at LambdaHat
 SclPrd = 0;
 for( Index i = 0 ; i < NumVar ; ++i )
  SclPrd += Gi[ i ] * ( Lambda[ i ] - LambdaHat[ i ] );

 SigmaHat += ( 1 - alpha ) * ( FiHat - FiLambda + SclPrd );

 // LambdaHat and its value
 comb( LambdaHat , 1 - alpha , Lambda , alpha );
 FiHat = alpha * FiLambda + ( 1 - alpha ) * FiHat;

 // the linearization error of the direction at LambdaHat
 HatEpsilon = alpha * SigmaHat + ( 1 - alpha ) * HatEpsilon;

 }  // end( SubGradSolver::chg_lambda_hat )

/*--------------------------------------------------------------------------*/

void SubGradSolver::update_cert_bound( void )
{
 // the aggregate linearization a^T x + c is a minorant of the function, and
 // its minimum on the box is a valid lower bound if it is finite
 const auto & a = deflection ? dirRaw : GiRaw;
 VarValue v = deflection ? dirCnst : GiCnst;

 for( Index i = 0 ; i < NumVar ; ++i )
  if( a[ i ] > 0 ) {
   if( lb[ i ] == - INF )
    return;
   v += a[ i ] * lb[ i ];
   }
  else
   if( a[ i ] < 0 ) {
    if( ub[ i ] == INF )
     return;
    v += a[ i ] * ub[ i ];
    }

 if( v > CertLB )
  CertLB = v;

 }  // end( SubGradSolver::update_cert_bound )

/*--------------------------------------------------------------------------*/

bool SubGradSolver::is_optimal( void ) const
{
 if( FiBest == INF )
  return( false );

 const VarValue eps = ( SGPar4 && ( ! NItIncr ) )
                      ? std::min( Epsilon , HatEpsilon ) : Epsilon;

 return( tStar * std::sqrt( NrmPrj ) + eps <= max_error( FiBest ) );
 }

/*--------------------------------------------------------------------------*/

VarValue SubGradSolver::lipschitz( void ) const
{
 VarValue L = 0;
 for( auto f : v_c05f )
  if( ( L += f->get_Lipschitz_constant() ) == INF )
   break;

 if( f_zero && ( L < INF ) )
  L += f_zero->get_Lipschitz_constant();

 // an infinite constant is estimated by the largest subgradient seen
 if( L == INF )
  L = std::sqrt( MaxNrmGi );

 return( L > 0 ? L : 1 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End CLASS SubGradSolver -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS Stepsize -------------------------------*/
/*--------------------------------------------------------------------------*/

Stepsize::StepsizeFactoryMap & Stepsize::f_factory( void )
{
 static StepsizeFactoryMap s_factory;
 return( s_factory );
 }

/*--------------------------------------------------------------------------*/

Stepsize::FunctionValue Stepsize::get_stepsize( bool incremental )
{
 static constexpr FunctionValue STEpZro = 1e-8;

 FunctionValue denomin;
 if( incremental ) {
  const auto L = read_Lipschitz();
  denomin = L * L * read_dbl_par( SubGradSolver::dblLipFactor ) *
            FunctionValue( read_NItIncr() );
  }
 else
  denomin = read_dir_norm();

 if( denomin <= STEpZro * std::max( std::abs( read_Fi_bar() ) ,
				    FunctionValue( 1 ) ) )
  return( 1 );

 auto Beta_ = Beta;
 if( Beta_ >= std::min( FunctionValue( 2 ) , MaxBeta ) )
  Beta_ = std::min( 2 - STEpZro , MaxBeta );

 return( Beta_ * ( read_Fi_bar() - FiLev ) / denomin );
 }

/*--------------------------------------------------------------------------*/

Stepsize::FunctionValue Stepsize::read_alpha( void ) const {
 return( Solver->alpha );
 }

Stepsize::FunctionValue Stepsize::read_Gi_norm( void ) const {
 return( Solver->NrmGi );
 }

Stepsize::FunctionValue Stepsize::read_dir_norm( void ) const {
 return( Solver->NrmDir );
 }

Stepsize::FunctionValue Stepsize::read_dGk( void ) const {
 return( Solver->dGk );
 }

Stepsize::FunctionValue Stepsize::read_dM1Gk( void ) const {
 return( Solver->dM1Gk );
 }

Stepsize::Index Stepsize::read_NItIncr( void ) const {
 return( Solver->NItIncr );
 }

Stepsize::FunctionValue Stepsize::read_Fi_lambda( void ) const {
 return( Solver->FiLambda );
 }

Stepsize::FunctionValue Stepsize::read_Fi_bar( void ) const {
 return( Solver->FiBar );
 }

Stepsize::FunctionValue Stepsize::read_Fi_best( void ) const {
 return( Solver->FiBest );
 }

Stepsize::FunctionValue Stepsize::read_lower_bound( void ) const {
 return( std::max( Solver->global_lb() , Solver->CondLB ) );
 }

Stepsize::FunctionValue Stepsize::read_Lipschitz( void ) const {
 return( Solver->lipschitz() );
 }

Stepsize::FunctionValue Stepsize::read_tStar( void ) const {
 return( Solver->tStar );
 }

int Stepsize::read_int_par( idx_type par ) const {
 return( Solver->get_int_par( par ) );
 }

double Stepsize::read_dbl_par( idx_type par ) const {
 return( Solver->get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Deflection ------------------------------*/
/*--------------------------------------------------------------------------*/

Deflection::DeflectionFactoryMap & Deflection::f_factory( void )
{
 static DeflectionFactoryMap s_factory;
 return( s_factory );
 }

/*--------------------------------------------------------------------------*/

Deflection::FunctionValue Deflection::read_stepsize( void ) const {
 return( Solver->step );
 }

Deflection::FunctionValue Deflection::read_Gi_norm( void ) const {
 return( Solver->NrmGi );
 }

Deflection::FunctionValue Deflection::read_dir_norm( void ) const {
 return( Solver->NrmDir );
 }

Deflection::FunctionValue Deflection::read_dGk( void ) const {
 return( Solver->dM1Gk );
 }

Deflection::FunctionValue Deflection::read_sigma( void ) const {
 return( Solver->Sigma );
 }

Deflection::FunctionValue Deflection::read_epsilon( void ) const {
 return( Solver->Epsilon );
 }

Deflection::FunctionValue Deflection::read_Fi_lambda( void ) const {
 return( Solver->FiLambda );
 }

Deflection::FunctionValue Deflection::read_Fi_bar( void ) const {
 return( Solver->FiBar );
 }

Deflection::FunctionValue Deflection::read_Fi_best( void ) const {
 return( Solver->FiBest );
 }

Deflection::Index Deflection::read_iterations( void ) const {
 return( Solver->TotIter );
 }

Deflection::FunctionValue Deflection::read_Lipschitz( void ) const {
 return( Solver->lipschitz() );
 }

Deflection::FunctionValue Deflection::read_tStar( void ) const {
 return( Solver->tStar );
 }

Deflection::FunctionValue Deflection::read_rel_acc( void ) const {
 return( Solver->RelAcc );
 }

int Deflection::read_int_par( idx_type par ) const {
 return( Solver->get_int_par( par ) );
 }

double Deflection::read_dbl_par( idx_type par ) const {
 return( Solver->get_dbl_par( par ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- End File SubGradSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
