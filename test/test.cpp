/*--------------------------------------------------------------------------*/
/*------------------------------ File test.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Unit test of SubGradSolver. Each SubGradSolver of the BlockSolverConfig
 * SGBSCfg.txt, one for each of the rules and of the schemes, solves in turn
 * three small problems whose optimum is known:
 *
 * - P1, the minimization of the convex PolyhedralFunction
 *   \f$ \max\{ | x_1 - 1 | , | x_2 - 2 | \} \f$ on \f$ R^2 \f$, with the
 *   valid lower bound 0, whose optimal value is 0;
 *
 * - P2, the minimization of \f$ \max\{ x_1 , x_2 , - x_1 - x_2 - 3 \} \f$
 *   on \f$ x \geq 0 \f$, whose optimal value is 0 at the origin, on the
 *   boundary of the box;
 *
 * - P3, the maximization of the sum of the linear function
 *   \f$ 0.5 x_1 - 0.5 x_2 \f$ and of the two concave functions
 *   \f$ \min\{ x_1 , 2 - x_1 \} \f$ and \f$ \min\{ 2 x_2 , 3 - x_2 \} \f$,
 *   each in its own sub-Block and depending on one ColVariable only, which
 *   is the structure of a Lagrangian dual; the optimal value is 3. The
 *   coefficient of \f$ x_1 \f$ is then changed to -0.5, which makes it 2,
 *   and \f$ x_2 \f$ is then fixed to 0.5, which makes it 1.25, so that the
 *   Modification are exercised.
 *
 * For each SubGradSolver and problem, the value of the point found must be
 * within the relative tolerance 1e-3 of the optimum if the status is kOK,
 * and 1e-2 if it is any other that is not an error (as a subgradient method
 * may well stop far from the optimum), the bounds must be valid, and in P3
 * the dual solution must be available. A table is printed with the values
 * and the statuses.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Matematica ed Informatica \n
 *         Universita' di Cagliari \n
 *
 * \author Laura Fiorella \n
 *         Dipartimento di Matematica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Gorgone, Laura Fiorella
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iomanip>
#include <iostream>

#include "AbstractBlock.h"
#include "BlockSolverConfig.h"
#include "FRealObjective.h"
#include "LinearFunction.h"
#include "PolyhedralFunction.h"

#include "SubGradSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using MultiVector = PolyhedralFunction::MultiVector;
using RealVector = PolyhedralFunction::RealVector;

/*--------------------------------------------------------------------------*/
/*------------------------------ CONSTANTS ---------------------------------*/
/*--------------------------------------------------------------------------*/

// the relative tolerance on the value of a SubGradSolver that claims to have
// found an optimal solution, and of one that has stopped for any other reason
static constexpr double TolOK = 1e-3;
static constexpr double TolStop = 1e-2;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/
// solve block with each of its SubGradSolver and compare with the optimum
// opt; dual says whether the dual solution has to be checked

static bool solve_all( Block * block , const std::string & name ,
		       double opt , bool max , bool dual )
{
 bool ok = true;
 int i = 0;
 for( auto s : block->get_registered_solvers() ) {
  ++i;
  auto sgs = dynamic_cast< SubGradSolver * >( s );
  if( ! sgs ) {
   std::cout << name << " #" << i << ": not a SubGradSolver" << std::endl;
   ok = false;
   continue;
   }

  const int st = sgs->compute( false );
  const double lb = sgs->get_lb();
  const double ub = sgs->get_ub();
  const double v = max ? lb : ub;
  const double tol = ( st == Solver::kOK ? TolOK : TolStop ) *
                     std::max( 1.0 , std::abs( opt ) );

  bool oki = ( ( ( st >= Solver::kOK ) && ( st < Solver::kError ) ) ||
	       ( st == Solver::kLowPrecision ) );
  oki = oki && ( std::abs( v - opt ) <= tol );
  oki = oki && ( lb <= opt + tol ) && ( ub >= opt - tol );

  if( dual && oki ) {
   oki = sgs->has_dual_solution();
   if( oki ) {
    sgs->get_dual_solution();
    for( SubGradSolver::Index k = 0 ; k < sgs->n_components() ; ++k ) {
     const auto & lc =
      sgs->component( k )->get_important_linearization_coefficients();
     if( ( lc.size() != 1 ) ||
	 ( ! sgs->component( k )->is_linearization_there( lc[ 0 ].first ) ) )
      oki = false;
     }
    }
   }

  std::cout << name << " #" << i << " " << std::setw( 10 )
	    << sgs->get_Stepsize()->classname() << " + " << std::setw( 10 )
	    << ( sgs->get_Deflection() ? sgs->get_Deflection()->classname()
		 : std::string( "none" ) )
	    << ": status " << std::setw( 2 ) << st << ", "
	    << std::setw( 6 ) << sgs->get_elapsed_iterations() << " it, value "
	    << std::scientific << std::setprecision( 6 ) << v << " (opt "
	    << opt << "), bounds [" << lb << " , " << ub << "] "
	    << ( oki ? "OK" : "KO" ) << std::endl;

  ok = ok && oki;
  }

 return( ok );
 }

/*--------------------------------------------------------------------------*/
// register the SubGradSolver of SGBSCfg.txt to block

static BlockSolverConfig * configure( Block * block )
{
 auto c = Configuration::deserialize( "SGBSCfg.txt" );
 auto bsc = dynamic_cast< BlockSolverConfig * >( c );
 if( ! bsc ) {
  delete c;
  throw( std::logic_error( "SGBSCfg.txt is not a BlockSolverConfig" ) );
  }

 bsc->apply( block );
 return( bsc );
 }

/*--------------------------------------------------------------------------*/
// unregister the SubGradSolver of bsc from block

static void unconfigure( Block * block , BlockSolverConfig * bsc )
{
 bsc->clear();
 bsc->apply( block );
 delete bsc;
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 bool ok = true;

 try {
  // P1: max{ | x_1 - 1 | , | x_2 - 2 | } on R^2 - - - - - - - - - - - - - -
  {
   std::vector< ColVariable > x( 2 );
   auto block = new AbstractBlock();
   block->add_static_variable( x , "x" );
   auto pf = new PolyhedralFunction( { & x[ 0 ] , & x[ 1 ] } ,
				     MultiVector( { { 1 , 0 } , { -1 , 0 } ,
						    { 0 , 1 } , { 0 , -1 } } ) ,
				     RealVector( { -1 , 1 , -2 , 2 } ) , 0 ,
				     true );
   auto obj = new FRealObjective( block , pf );
   obj->set_sense( Objective::eMin , eNoMod );
   block->set_objective( obj , eNoMod );

   auto bsc = configure( block );
   ok = solve_all( block , "P1" , 0 , false , false ) && ok;
   unconfigure( block , bsc );
   delete block;  // which does not delete the Objective
   delete obj;
   }

  // P2: max{ x_1 , x_2 , - x_1 - x_2 - 3 } on x >= 0- - - - - - - - - - - -
  {
   std::vector< ColVariable > x( 2 );
   for( auto & xi : x ) {
    xi.is_positive( true , eNoMod );
    xi.set_value( 5 );
    }
   auto block = new AbstractBlock();
   block->add_static_variable( x , "x" );
   auto pf = new PolyhedralFunction( { & x[ 0 ] , & x[ 1 ] } ,
				     MultiVector( { { 1 , 0 } , { 0 , 1 } ,
						    { -1 , -1 } } ) ,
				     RealVector( { 0 , 0 , -3 } ) , -1 ,
				     true );
   auto obj = new FRealObjective( block , pf );
   obj->set_sense( Objective::eMin , eNoMod );
   block->set_objective( obj , eNoMod );

   auto bsc = configure( block );
   ok = solve_all( block , "P2" , 0 , false , false ) && ok;
   unconfigure( block , bsc );
   delete block;
   delete obj;
   }

  // P3: a Lagrangian-like concave maximization with sub-Block- - - - - - - -
  {
   std::vector< ColVariable > x( 2 );
   auto root = new AbstractBlock();
   root->add_static_variable( x , "x" );

   auto lf = new LinearFunction( LinearFunction::v_coeff_pair(
				   { { & x[ 0 ] , 0.5 } , { & x[ 1 ] , -0.5 } }
				   ) );
   auto robj = new FRealObjective( root , lf );
   robj->set_sense( Objective::eMax , eNoMod );
   root->set_objective( robj , eNoMod );

   const MultiVector A[ 2 ] = { { { 1 } , { -1 } } , { { 2 } , { -1 } } };
   const RealVector b[ 2 ] = { { 0 , 2 } , { 0 , 3 } };
   std::vector< FRealObjective * > sbobj;
   for( int k = 0 ; k < 2 ; ++k ) {
    auto sb = new AbstractBlock( root );
    auto pf = new PolyhedralFunction( { & x[ k ] } , MultiVector( A[ k ] ) ,
				      RealVector( b[ k ] ) ,
				      Inf< double >() , false );
    auto obj = new FRealObjective( sb , pf );
    obj->set_sense( Objective::eMax , eNoMod );
    sb->set_objective( obj , eNoMod );
    root->add_nested_Block( sb );
    sbobj.push_back( obj );
    }

   // a valid upper bound, which Polyak and FumeroTV need
   root->set_valid_upper_bound( 3.5 );

   auto bsc = configure( root );
   ok = solve_all( root , "P3" , 3 , true , true ) && ok;

   // change the coefficient of x_1: the optimum becomes 2
   lf->modify_coefficient( 0 , -0.5 );
   ok = solve_all( root , "P3'" , 2 , true , true ) && ok;

   // fix x_2 to 0.5: the optimum becomes 1.25
   x[ 1 ].set_value( 0.5 );
   x[ 1 ].is_fixed( true );
   ok = solve_all( root , "P3''" , 1.25 , true , true ) && ok;

   unconfigure( root , bsc );
   delete root;
   delete robj;
   for( auto obj : sbobj )
    delete obj;
   }
  }
 catch( std::exception & e ) {
  std::cerr << "SubGradSolver_unit_test: " << e.what() << std::endl;
  return( 1 );
  }

 std::cout << "SubGradSolver: " << ( ok ? "all tests passed" :
				     "some tests FAILED" ) << std::endl;

 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- End File test.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
