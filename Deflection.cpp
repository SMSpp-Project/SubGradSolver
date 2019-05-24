/*
 * Deflection.cpp
 *
 *  Created on: May 10, 2019
 *      Author: simone
 */




#include "Deflection.h"


SMSpp_insert_in_factory_cpp_1(PrimalDual);
SMSpp_insert_in_factory_cpp_1(Volume);



/*--------------------------------------------------------------------------*/
/*------------------------- METHODS of Deflection --------------------------*/
/*--------------------------------------------------------------------------*/

Deflection::Deflection(SUBGRADSolver* slvr) {

	if( ! slvr )
	     throw(std::logic_error( "Deflection::Deflection: no subgradient solver " ));
	    Solver = slvr;

	    VOLLog = 0;
	    VOLLLvl = 0;


	    auto extraconf=slvr->get_ComputeConfig(true);

	    SimpleConfiguration<std::pair< Configuration * , Configuration * >> *sc =dynamic_cast<SimpleConfiguration<std::pair< Configuration * , Configuration * >>*> (extraconf->f_extra_Configuration);

	    if(sc==nullptr)
	    	throw( std::logic_error( "the extra configuration is not a SimpleConfiguration<std::pair< Configuration * , Configuration * >> " ) );

		deflectionconf = dynamic_cast<SimpleConfiguration<std::vector<double>> *> (sc->f_value.second);

		if(deflectionconf==nullptr)
			throw( std::logic_error( "the  configuration of deflection is not a SimpleConfiguration<std::vector<double>> " ) );




}

Deflection* Deflection::new_Deflection(const std::string& classname,
		SUBGRADSolver* slvr) {

	const auto it = Deflection::f_factory().find( classname );
			   if( it == Deflection::f_factory().end() )
			    throw( std::invalid_argument( classname +
			 			   std::string( " not present in Stepsize factory" ) ) );

			   return ((it->second)(slvr));
}



Deflection::DeflectionFactoryMap& Deflection::f_factory(void) {

	static DeflectionFactoryMap s_factory;
	 return( s_factory );
}


void Deflection::SetVOLLog(std::ostream* outs, const char lvl) {

	if( ( VOLLog = outs ) )
	     VOLLLvl = lvl;
	    else
	     VOLLLvl = 0;

}


double Deflection::Delta(void) {
	return( Inf<double>() );

}


double Deflection::GetStepsize(void) {
	return this->Solver->step;
}

const double Deflection::GetGiNorm(void) {
	return this->Solver->NrmGi;
}

const double Deflection::GetDNorm(void) {
	return this->Solver->NrmDir;
}

const double Deflection::GetdGk(void) {
	return this->Solver->dGk;
}

const double Deflection::GetSigma(void) {
	return this->Solver->Sigma;
}

const double Deflection::GetEpsilon(void) {
	return this->Solver->Epsilon;
}

void Deflection::Format(void) {
}

const bool Deflection::DoSS(void) {
	return( true );
}

const double Deflection::ReadFVal(void) {
	return this->Solver->FiLambda;
}





/*--------------------------------------------------------------------------*/
/*------------------------- METHODS of PrimalDual --------------------------*/
/*--------------------------------------------------------------------------*/

PrimalDual::PrimalDual(SUBGRADSolver* slvr) : Deflection::Deflection( slvr ),Stepsize::Stepsize(slvr){

	this->average=deflectionconf->f_value[par_average];


}

void PrimalDual::SetVOLLog(std::ostream* outs, const char lvl) {

	Deflection::SetVOLLog( outs , lvl );

	 #if( LOG_VOL )
	  if( VOLLLvl > 1 )
	   if( average )
	    *VOLLog << endl << "PrimalDual: Simple averages ";
	   else
	    *VOLLog << endl << "PrimalDual: Weighted averages ";

	  *VOLLog << "~ GammaFact = " <<  LpsFct << endl;
	 #endif


}

void PrimalDual::Format(void) {

	Deflection::Format();

	 double LipTmp = GetOracle()->get_Lipschitz_constant();
	 if( LipTmp != LipCnst )
	  LipCnst = LipTmp;

	 beta1 = 1;
	 delta = 0;
}

void PrimalDual::NewDEF(void) {
	if( average )
	  vk = 1;
	 else
	  vk = 1 / std::sqrt( Deflection::GetGiNorm() );

	 Tmpdelta = delta + vk;


}

void PrimalDual::NewStep(void) {
	// update the Lipschitz constant - - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 double LipTmp = GetOracle()->get_Lipschitz_constant();
	 if( LipTmp != LipCnst )
	  LipCnst = LipTmp;

	 UpdateGamma( );

	 // SR computation - - -  - - - - - - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 StepSize = Tmpdelta / ( gamma * beta1 );

	 VOLLOG( 1, endl << "            Delta = " << Tmpdelta << " ~ Beta = "
		    << beta1 << " ~ gamma = " << gamma );

	 delta = Tmpdelta;
	 beta1 = ( beta1 + 1/beta1 );
}

double PrimalDual::GetDFLCoeff(void) {
	return( vk / Tmpdelta );
}

double PrimalDual::GetStepsize(bool StepIsIncr) {
	if( StepIsIncr )
	  throw( std::logic_error( "PrimalDual::GetStepsize: this is not allowed" ) );

	 return( StepSize );
}

void PrimalDual::UpdateGamma(void) {
	 double tStar;

	 tStar = Deflection::Solver->get_dbl_par(Deflection::Solver->dbl_par_str2idx("ktStar"));


	 if( average )
	  gamma = LpsFct * ( std::sqrt( LipCnst / 2.0 )  / tStar );
	 else
	  gamma = LpsFct / ( std::sqrt( 2.0 * LipCnst ) * tStar );

}

double PrimalDual::GetLev(void) {
	throw( std::logic_error( "PrimalDual::GetLev: this call is not allowed" ) );
}

double PrimalDual::GetBeta(void) {
	throw( std::logic_error( "PrimalDual::GetBeta: this call is not allowed" ) );
}

void PrimalDual::SetMaxBeta(const double alpha) {
	throw( std::logic_error( "PrimalDual::SetMaxBeta: this call is not allowed" ) );
}



/*--------------------------------------------------------------------------*/
/*------------------------- METHODS of Volume ------------------------------*/
/*--------------------------------------------------------------------------*/

Volume::Volume(SUBGRADSolver* slvr) : Deflection::Deflection( slvr ){


	this->tauInit=deflectionconf->f_value[par_tauInit]; //
	this->tauMin=deflectionconf->f_value[par_tauMin];
	this->tauFactor=deflectionconf->f_value[par_tauFactor];
	this->tauIter=deflectionconf->f_value[par_tauIter];
	this->m=deflectionconf->f_value[par_m];

	alpha = lastvalue=tau=double(0);

}

void Volume::SetVOLLog(std::ostream* outs, const char lvl) {

	Deflection::SetVOLLog( outs , lvl );

	 #if( LOG_VOL )
	  if( VOLLLvl > 1 )
	   *VOLLog << std::endl <<  "Volume: tau = " << tauInit << " in ["<< tauMin
		   << " , -] with ( " << tauFactor << " , " << tauIter << " ) ~ m = "
		   << m << std::endl;
	 #endif

}

void Volume::Format(void) {
	Deflection::Format( );
	alpha = lastvalue = Inf<double>();  // deflection coefficient is unknown
	tau = tauInit;
}

void Volume::NewDEF(void) {

	// change tau if needed - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( Solver->NrIter() <= 1 ) {
	  alpha = 1;
	  return;
	  }

	 double FiRec = Solver->ReadBestFiVal();

	 // tau should be decreased if after tauIter iterations Fi() is decreased
	 // less than 1% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( ! ( Solver->NrIter() % tauIter ) )
	  if( tau > tauMin ) {
	   const double x = ( lastvalue - FiRec ) / ABS( FiRec );
	   lastvalue = FiRec;
	   if( x <= 0.01 ) {
	    tau *= tauFactor;
	    if( tau < tauMin )
	     tau = tauMin;
	    }
	   }

	 // solve a quadratic problem and get the deflection coefficient - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 Solve();

}

double Volume::GetDFLCoeff(void) {
	return( alpha );
}

const bool Volume::DoSS(void) {

	bool serious;
	double Delta;

	Delta = Solver->get_dbl_par(Solver->dbl_par_str2idx("par_EpsLin"));
	Delta *= m * std::max( 1.0 , ABS( Solver->ReadBestFiVal() ) );

	 // if there is an improvement of Fi(), a SS occurs; otherwise, a NS occurs
	 if( ( Solver->ReadFiVal() - Solver->ReadFVal() ) >= Delta )
	  serious = true;
	 else
	  serious = false;

	 return( serious );
}

double Volume::Delta(void) {
	double EpsLin;

	EpsLin = Solver->get_dbl_par(Solver->dbl_par_str2idx("par_EpsLin"));

	 return( m * EpsLin * std::max( 1.0 , ABS( Solver->ReadBestFiVal() ) ) );
}

void Volume::Solve(void) {

	double GiNorm = GetGiNorm();
	double DNorm = GetDNorm();
	double GixD = GetdGk();
	double Sigma = GetSigma();
	double Epsilon = GetEpsilon();

	 // compute the proximity parameter  - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 double proximity = GetStepsize() == Inf<double>() ? 1 : GetStepsize();

	 // solving the cutting plane model  - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 double den = proximity * ( GiNorm + DNorm - 2.0 * GixD );

	 if( den <= 1e-8 )
	  alpha = 1;
	 else
	  alpha = ( Epsilon - Sigma - proximity * ( GixD - DNorm ) ) / den;

	 // tau should be projected over [0,1] - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( alpha > 1 )
	  alpha = std::min( tau , 1.0 );
	 else
	  if( alpha < 1e-8 )
	   alpha = tau / 10;

}


