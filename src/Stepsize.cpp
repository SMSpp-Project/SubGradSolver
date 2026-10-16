/*--------------------------------------------------------------------------*/
/*--------------------------- File Stepsize.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the stepsize rules Polyak, ColorTV and FumeroTV [see
 * Polyak.h, ColorTV.h and FumeroTV.h], and their registration in the
 * Stepsize factory.
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ColorTV.h"

#include "FumeroTV.h"

#include "Polyak.h"

#include "SubGradSolver.h"

#include <cmath>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( Polyak );

SMSpp_insert_in_factory_cpp_0( ColorTV );

SMSpp_insert_in_factory_cpp_0( FumeroTV );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF Polyak ----------------------------*/
/*--------------------------------------------------------------------------*/

void Polyak::format( void )
{
 Stepsize::format();
 Beta = read_dbl_par( SubGradSolver::dblPolyakBeta );
 }

/*--------------------------------------------------------------------------*/

void Polyak::new_step( void )
{
 Beta = read_dbl_par( SubGradSolver::dblPolyakBeta );

 update_target_level();

 if( FiLev == - Inf< FunctionValue >() )
  throw( std::logic_error( "Polyak::new_step: no lower bound is known" ) );
 }

/*--------------------------------------------------------------------------*/

bool Polyak::update_target_level( void )
{
 const auto LwrBnd = read_lower_bound();
 if( LwrBnd > FiLev ) {
  FiLev = LwrBnd;
  return( true );
  }

 return( false );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- METHODS OF ColorTV ----------------------------*/
/*--------------------------------------------------------------------------*/

void ColorTV::format( void )
{
 Stepsize::format();
 lastgreeniter = lastyellowiter = lastrediter = 0;
 LwrBnd = - Inf< FunctionValue >();
 lastvalue = FiLmbd = Inf< FunctionValue >();
 Beta = read_dbl_par( SubGradSolver::dblColorBeta0 );
 NrIter = 0;
 }

/*--------------------------------------------------------------------------*/

void ColorTV::new_step( void )
{
 FiLmbd = read_Fi_lambda();
 if( FiLmbd == Inf< FunctionValue >() )
  throw( std::logic_error( "ColorTV::new_step: unknown function value" ) );

 update_target_level();

 coloring();

 ++NrIter;
 }

/*--------------------------------------------------------------------------*/

bool ColorTV::update_target_level( void )
{
 bool TLchgd = false;

 const auto LwrBnd_ = read_lower_bound();
 if( LwrBnd_ > LwrBnd )
  LwrBnd = LwrBnd_;

 if( FiLev == - Inf< FunctionValue >() ) {
  // no level yet: start 5% below the current value
  FiLev = FiLmbd - 0.05 * std::max( std::abs( FiLmbd ) ,
				    FunctionValue( 1 ) );
  TLchgd = true;
  }
 else
  if( ( FiLev < -10 ) || ( FiLev > 10 ) ) {  // nicely away from 0
   if( FiLmbd <= FiLev + std::abs( FiLev ) * 0.05 ) {
    if( ( FiLmbd < 10 ) && ( FiLmbd > -9.5 ) )
     FiLev = -10;
    else {
     FiLev -= 0.025 * std::abs( FiLev );
     FiLev = std::min( FiLev , FiLmbd - 0.05 * std::abs( FiLmbd ) );
     }
    TLchgd = true;
    }
   }
  else                                       // near 0, must be careful
   if( FiLmbd < 10 ) {
    FiLev = ( FiLmbd < -9.5 ) ? -10 :
            ( FiLmbd > 0 ? FiLmbd * 0.95 : FiLmbd * 1.05 );
    TLchgd = true;
    }

 // the level is never below the lower bound
 if( FiLev < LwrBnd ) {
  FiLev = LwrBnd;
  TLchgd = true;
  }

 return( TLchgd );

 }  // end( ColorTV::update_target_level )

/*--------------------------------------------------------------------------*/

ColorTV::condition ColorTV::coloring( void )
{
 const Index greentestinvl = read_int_par( SubGradSolver::intColorGreen );
 const Index yellowtestinvl = read_int_par( SubGradSolver::intColorYellow );
 const Index redtestinvl = read_int_par( SubGradSolver::intColorRed );

 condition lastswing;
 if( ( read_dGk() >= 1e-6 ) &&
     ( lastvalue - FiLmbd >= 1e-6 * std::max( std::abs( read_Fi_best() ) ,
					     FunctionValue( 1 ) ) ) ) {
  lastswing = green;
  lastgreeniter = NrIter;
  }
 else
  if( ( read_dGk() < 1e-6 ) && ( lastvalue > FiLmbd ) ) {
   lastswing = yellow;
   lastyellowiter = NrIter;
   }
  else {
   lastswing = red;
   lastrediter = NrIter;
   }

 // change beta, if needed
 switch( lastswing ) {
  case( green ):
   if( ( NrIter - std::max( lastyellowiter , lastrediter ) >=
	 greentestinvl ) && ( Beta < 2 ) ) {
    lastgreeniter = lastyellowiter = lastrediter = NrIter;
    Beta *= 2;
    }
   break;
  case( yellow ):
   if( NrIter - std::max( lastgreeniter , lastrediter ) >= yellowtestinvl ) {
    lastgreeniter = lastyellowiter = lastrediter = NrIter;
    Beta *= 1.1;
    }
   break;
  case( red ):
   if( ( NrIter - std::max( lastgreeniter , lastyellowiter ) >=
	 redtestinvl ) && ( Beta > FunctionValue( 5e-4 ) ) ) {
    lastgreeniter = lastyellowiter = lastrediter = NrIter;
    Beta *= 0.67;
    }
   break;
  }

 lastvalue = read_Fi_bar();

 return( lastswing );

 }  // end( ColorTV::coloring )

/*--------------------------------------------------------------------------*/
/*------------------------- METHODS OF FumeroTV ----------------------------*/
/*--------------------------------------------------------------------------*/

void FumeroTV::format( void )
{
 Stepsize::format();
 Beta = read_dbl_par( SubGradSolver::dblFumeroBeta0 );
 FiLambda = FiRef = Inf< FunctionValue >();
 LwrBnd = - Inf< FunctionValue >();
 etac = rc = 0;
 const auto sigmaMin = read_dbl_par( SubGradSolver::dblFumeroSigmaMin );
 const auto r1 = read_int_par( SubGradSolver::intFumeroR1 );
 r2 = Index( std::ceil( std::pow( - std::log( sigmaMin ) / 0.6933 ,
				  1.0 / 3.26 ) * r1 ) );
 sigma = 1;
 }

/*--------------------------------------------------------------------------*/

void FumeroTV::new_step( void )
{
 FiLambda = read_Fi_lambda();
 if( FiLambda == Inf< FunctionValue >() )
  throw( std::logic_error( "FumeroTV::new_step: unknown function value" ) );

 if( update_target_level() )
  FiLev = sigma * LwrBnd + ( 1 - sigma ) * FiRef;

 if( Beta < 1e-6 )
  Beta = 1e-6;
 }

/*--------------------------------------------------------------------------*/

bool FumeroTV::update_target_level( void )
{
 bool TLchgd = false;

 const auto LwrBnd_ = read_lower_bound();
 if( LwrBnd_ > LwrBnd ) {
  LwrBnd = LwrBnd_;
  TLchgd = true;
  }

 if( LwrBnd == - Inf< FunctionValue >() )
  throw( std::logic_error(
		  "FumeroTV::update_target_level: no lower bound is known" ) );

 const Index eta1 = read_int_par( SubGradSolver::intFumeroEta1 );
 const Index eta2 = read_int_par( SubGradSolver::intFumeroEta2 );

 const auto FiRef_ = ( FiRef == Inf< FunctionValue >() ) ? FiRef :
  FiRef - 1e-6 * std::max( FunctionValue( 1 ) , std::abs( FiRef ) );

 if( FiLambda <= FiRef_ ) {  // an improvement
  FiRef = FiLambda;
  if( rc >= r2 )
   Beta *= 2;
  etac = 0;
  TLchgd = true;
  }
 else                        // a failure
  if( rc >= r2 ) {           // phase two
   if( ++etac >= eta1 ) {
    Beta /= 2;
    etac = 0;
    }
   }
  else                       // phase one
   if( ++etac >= eta2 ) {
    Beta = Beta / ( 2 * Beta + 1 );
    sigma = exp_fun( ++rc );
    TLchgd = true;
    etac = 0;
    }

 if( FiLev < LwrBnd ) {
  FiLev = LwrBnd;
  TLchgd = true;
  }

 return( TLchgd );

 }  // end( FumeroTV::update_target_level )

/*--------------------------------------------------------------------------*/

FumeroTV::FunctionValue FumeroTV::exp_fun( Index ri ) const
{
 const auto r1 = read_int_par( SubGradSolver::intFumeroR1 );
 return( std::exp( -0.6933 * std::pow( FunctionValue( ri ) /
				       FunctionValue( r1 ) , 3.26 ) ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File Stepsize.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
