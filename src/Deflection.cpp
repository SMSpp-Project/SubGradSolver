/*--------------------------------------------------------------------------*/
/*-------------------------- File Deflection.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the deflection rules Volume and PrimalDual [see
 * Volume.h and PrimalDual.h], and their registration in the Deflection
 * factory (and, for PrimalDual, also in the Stepsize one).
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

#include "PrimalDual.h"

#include "SubGradSolver.h"

#include "Volume.h"

#include <cmath>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( Volume );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// PrimalDual is registered in both factories, which the macros cannot do

SMSpp_define_force_load( PrimalDual )

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
 }

PrimalDual::_init PrimalDual::_initializer;

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF Volume ----------------------------*/
/*--------------------------------------------------------------------------*/

void Volume::format( void )
{
 Deflection::format();
 alpha = 1;
 lastvalue = Inf< FunctionValue >();
 tau = read_dbl_par( SubGradSolver::dblVolumeTauInit );
 }

/*--------------------------------------------------------------------------*/

void Volume::new_deflection( void )
{
 const auto FiRec = read_Fi_best();
 const Index tauIter = read_int_par( SubGradSolver::intVolumeTauIter );
 const auto tauMin = read_dbl_par( SubGradSolver::dblVolumeTauMin );

 // tau is decreased if in the last tauIter iterations the best value has
 // decreased by less than 1%
 if( tauIter && read_iterations() && ( ! ( read_iterations() % tauIter ) ) &&
     ( tau > tauMin ) ) {
  if( ( lastvalue < Inf< FunctionValue >() ) &&
      ( ( lastvalue - FiRec ) / std::max( std::abs( FiRec ) ,
					  FunctionValue( 1e-12 ) ) <= 0.01 ) )
   tau = std::max( tau * read_dbl_par( SubGradSolver::dblVolumeTauFactor ) ,
		   tauMin );
  lastvalue = FiRec;
  }

 // solve the one-dimensional quadratic problem
 const auto GiNorm = read_Gi_norm();
 const auto DNorm = read_dir_norm();
 const auto GixD = read_dGk();
 const auto Sigma = read_sigma();
 const auto Epsilon = read_epsilon();

 const auto stp = read_stepsize();
 const auto proximity = ( ( stp <= 0 ) || ( stp == Inf< FunctionValue >() ) )
                        ? 1 : stp;

 const auto den = proximity * ( GiNorm + DNorm - 2 * GixD );
 if( den <= 1e-8 )
  alpha = 1;
 else
  alpha = ( Epsilon - Sigma - proximity * ( GixD - DNorm ) ) / den;

 // safeguard alpha in [ tau / 10 , tau ]
 if( alpha > 1 )
  alpha = std::min( tau , FunctionValue( 1 ) );
 else
  if( alpha < 1e-8 )
   alpha = tau / 10;

 }  // end( Volume::new_deflection )

/*--------------------------------------------------------------------------*/

bool Volume::do_SS( void )
{
 return( read_Fi_bar() - read_Fi_lambda() >= delta() );
 }

/*--------------------------------------------------------------------------*/

Volume::FunctionValue Volume::delta( void ) const
{
 return( read_dbl_par( SubGradSolver::dblVolumeM ) * read_rel_acc() *
	 std::max( FunctionValue( 1 ) , std::abs( read_Fi_best() ) ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- METHODS OF PrimalDual --------------------------*/
/*--------------------------------------------------------------------------*/

void PrimalDual::format( void )
{
 Deflection::format();
 Stepsize::format();
 LipCnst = Deflection::read_Lipschitz();
 beta1 = 1;
 delta_ = 0;
 tmp_delta = 1;
 vk = 1;
 }

/*--------------------------------------------------------------------------*/

void PrimalDual::new_deflection( void )
{
 if( Deflection::read_int_par( SubGradSolver::intPDAverage ) )
  vk = 1;
 else {
  const auto n = std::sqrt( Deflection::read_Gi_norm() );
  vk = ( n > 0 ) ? 1 / n : 1;
  }

 tmp_delta = delta_ + vk;
 }

/*--------------------------------------------------------------------------*/

void PrimalDual::new_step( void )
{
 LipCnst = Deflection::read_Lipschitz();

 update_gamma();

 step_ = tmp_delta / ( gamma * beta1 );

 delta_ = tmp_delta;
 beta1 = beta1 + 1 / beta1;
 }

/*--------------------------------------------------------------------------*/

PrimalDual::FunctionValue PrimalDual::get_stepsize( bool incremental )
{
 if( incremental )
  throw( std::logic_error( "PrimalDual::get_stepsize: the incremental "
			   "variant is not supported" ) );
 return( step_ );
 }

/*--------------------------------------------------------------------------*/

PrimalDual::FunctionValue PrimalDual::get_level( void ) const
{
 throw( std::logic_error( "PrimalDual::get_level: there is no level" ) );
 }

/*--------------------------------------------------------------------------*/

PrimalDual::FunctionValue PrimalDual::get_beta( void ) const
{
 throw( std::logic_error( "PrimalDual::get_beta: there is no beta" ) );
 }

/*--------------------------------------------------------------------------*/

void PrimalDual::set_max_beta( FunctionValue alpha )
{
 throw( std::logic_error( "PrimalDual::set_max_beta: there is no beta" ) );
 }

/*--------------------------------------------------------------------------*/

void PrimalDual::update_gamma( void )
{
 const auto tStar = Deflection::read_tStar();
 const auto F = Deflection::read_dbl_par( SubGradSolver::dblLipFactor );

 if( Deflection::read_int_par( SubGradSolver::intPDAverage ) )
  gamma = F * ( std::sqrt( LipCnst / 2 ) / tStar );
 else
  gamma = F / ( std::sqrt( 2 * LipCnst ) * tStar );

 if( ! ( gamma > 0 ) )  // a degenerate factor: no step is better than INF
  gamma = Inf< FunctionValue >();
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File Deflection.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
