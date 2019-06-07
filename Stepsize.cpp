/*
 * Stepsize.cpp
 *
 *  Created on: May 4, 2019
 *      Author: simone
 */


#include "Stepsize.h"





SMSpp_insert_in_factory_cpp_1(ColorTV);
SMSpp_insert_in_factory_cpp_1(FumeroTV);
SMSpp_insert_in_factory_cpp_1(Polyak);



/*--------------------------------------------------------------------------*/
/*------------------------- METHODS of Stepsize ----------------------------*/
/*--------------------------------------------------------------------------*/



Stepsize* Stepsize::new_Stepsize(const std::string& classname,
		SubGradSolver* slvr) {

	const auto it = Stepsize::f_factory().find( classname );
		   if( it == Stepsize::f_factory().end() )
		    throw( std::invalid_argument( classname +
		 			   std::string( " not present in Stepsize factory" ) ) );

		   return ((it->second)(slvr));


}

Stepsize::StepsizeFactoryMap& Stepsize::f_factory(void) {

	static StepsizeFactoryMap s_factory;
	 return( s_factory );
}




Stepsize::Stepsize(SubGradSolver* slvr) {

	Solver=slvr; ///< (pointer to) the SubGrad solver


	STPLog=nullptr;  ///< the output stream object
	STPLLvl = char(0);   ///< the "level of verbosity"

	FiLev=double(0);      ///< the target level \f$ f^{lev}_i \f$
	Beta= double(0);     ///< beta factor \f$ \beta_i \f$
	MaxBeta=double(0);     ///< maximum value for beta factor

	auto extraconf=slvr->get_ComputeConfig(true);

	SimpleConfiguration<std::pair< Configuration * , Configuration * >> *sc =dynamic_cast<SimpleConfiguration<std::pair< Configuration * , Configuration * >>*> (extraconf->f_extra_Configuration);

	if(sc==nullptr)
		 throw( std::logic_error( "the extra configuration is not a SimpleConfiguration<std::pair< Configuration * , Configuration * >> " ) );

	stepsizeconf = dynamic_cast<SimpleConfiguration<std::vector<double>> *> (sc->f_value.first);

	if(stepsizeconf==nullptr)
		throw( std::logic_error( "the  configuration of stepsize is not a SimpleConfiguration<std::vector<double>> " ) );




	 LpsFct=stepsizeconf->f_value[par_LpsFct];  ///< scaling factor



}

void Stepsize::SetSTPLog(std::ostream* outs, const char lvl) {

	if( ( STPLog = outs ) )
	     STPLLvl = lvl;
	    else
	     STPLLvl = 0;
}

void Stepsize::Format(void) {
	Beta = MaxBeta = Inf<double>();
	FiLev = - Inf<double>();
}

bool Stepsize::NeedsdkM1Gk(void) {
	return( false );
}

double Stepsize::GetStepsize(bool StepIsIncr) {

	double denomin, step;
	const double STEpZro = 1e-8;

	if( StepIsIncr ) {

	 Function::FunctionValue  LipCns =0;
	 for(auto a: GetOracle()){
		 LipCns+=a->get_Lipschitz_constant();
	 }

	 denomin = LipCns * LipCns * Function::FunctionValue( LpsFct * GetNItIcr() );
	 }
	else
	 denomin = GetDNorm();

	if( denomin <= ( STEpZro * std::max( std::abs( ReadFiBar() ) , Function::FunctionValue(1) ) ) )
	 step = 1;
	else {
		Function::FunctionValue Beta_ = Beta;
	 if( Beta_ >= std::min( 2.0 , MaxBeta ) )
	  Beta_ = std::min( 2.0 - STEpZro , MaxBeta );
	 step = Beta_ * ( ReadFiBar() - FiLev ) / denomin;
	 }

	return( step );


}

double Stepsize::GetLev(void) {
	return( FiLev );
}

double Stepsize::GetBeta(void) {
	return Beta;
}

void Stepsize::SetMaxBeta(const double alpha) {
	 MaxBeta = alpha;
}

bool Stepsize::UpdateTargetLevel(void) {
	return false;
}





double Stepsize::GetCoeffDefl(void) {
	return this->Solver->alpha;
}

std::vector<C05Function*> Stepsize::GetOracle(void) {
	return this->Solver->v_c05f;
}

const double Stepsize::GetGiNorm(void) {
	return this->Solver->NrmGi;
}

const double Stepsize::GetDNorm(void) {
	return this->Solver->NrmDir;
}

const double Stepsize::GetdGk(void) {
	return this->Solver->dGk;
}

double Stepsize::GetdkM1Gk(void) {
	return this->Solver->dM1Gk;
}

Variable::Index Stepsize::GetNItIcr(void) {
	return this->Solver->NItIncr;
}

const double Stepsize::ReadFkVal(void) {
	return this->Solver->FiLambda;
}

const double Stepsize::ReadFiBar(void) {
	return this->Solver->FiBar;
}



/*--------------------------------------------------------------------------*/
/*------------------------- METHODS of ColorTV -----------------------*/
/*--------------------------------------------------------------------------*/


ColorTV::ColorTV(SubGradSolver* slvr) : Stepsize::Stepsize( slvr  ) {

	this->BetaZero=stepsizeconf->f_value[par_BetaZero]; // initial value of beta
	this->greentestinvl=stepsizeconf->f_value[par_greentestinvl];
	this->yellowtestinvl=stepsizeconf->f_value[par_yellowtestinvl];
	this->redtestinvl=stepsizeconf->f_value[par_redtestinvl];





}

void ColorTV::SetSTPLog(std::ostream* outs, const char lvl) {
	Stepsize::SetSTPLog( outs , lvl );

	 #if( LOG_STP )
	  if( STPLLvl > 1 ) {
	   *STPLog << std::endl << "ColorTV: ~ BetaZero = " << BetaZero
		   << " ~ Color (g , y , r) = (" << greentestinvl << " , "
		   << yellowtestinvl << " , " << redtestinvl << ")" << std::endl;
	   }
	 #endif

}

void ColorTV::Format(void) {

	Stepsize::Format();

	lastgreeniter = lastyellowiter = lastrediter = 0;
	LwrBnd = - Inf<double>();   // lower bound, target level are unknown
	lastvalue = FiLmbd = Inf<double>(); // \f$ \bar{f}_{i-1} \f$ and \f$ f_i \f$
	                                    // are unknown
	Beta = BetaZero;
	NrIter = 0;

}

void ColorTV::NewStep(void) {

	// get the function value - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 FiLmbd = ReadFkVal();

	 if( FiLmbd == Inf<double>() )
		 std::logic_error( "ColorTV::GetStepsize: this should not happen" );

	 // eventually, improve the target level  - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 UpdateTargetLevel();

	 // which is it the color?  - - - - - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 coloring();
	 NrIter++;

	 // print the level   - - - - - - - - - - - - - - - - - - - - - - - - - - - -


	 //log comment
//	 STPLOG( 1 , std::endl << "           " << " level  = " << -FiLev
//		               << " ~  beta = " << Beta << std::endl << "           "
//		 );

}

bool ColorTV::UpdateTargetLevel(void) {

	bool TLchgd = false;

	 // the target level is updated after a change of the lower bound  - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


	 Function::FunctionValue LwrBnd_ = 0;
	 for(auto a : GetOracle()){
		 LwrBnd_+=a->get_lower_estimate();
	 }


	 if( LwrBnd_ > LwrBnd )
	  LwrBnd = LwrBnd_;

	 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( LwrBnd == - Inf<double>() )
		 std::logic_error("ColorTV::UpdateTargetLevel: no lower bound is given" );

	 // target level updating  - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( ( FiLev < -10 ) || ( FiLev > 10 ) ) {  // nicely away from 0
	  if( FiLmbd <= FiLev + std::abs( FiLev ) * 0.05 ) {
	   if( ( FiLmbd < 10 ) && ( FiLmbd > - 9.5 ) )
	    FiLev= - 10;
	   else {
	    FiLev -= 0.025 * std::abs( FiLev  );
	    FiLev = std::min( FiLev , FiLmbd - 0.05 * std::abs( FiLmbd ) );
	    }
	   TLchgd = true;
	   }
	  }
	 else                                       // near 0, must be careful
	  if( FiLmbd < 10 ) {
	   FiLev = ( FiLmbd < - 9.5 ) ? - 10 :
	            ( FiLmbd > 0 ? FiLmbd * 0.95 : FiLmbd * 1.05 );
	   TLchgd = true;
	   }

	 // failure of the target level updating - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( FiLev < LwrBnd ) {  // change the level
	  FiLev = LwrBnd;
	  TLchgd = true;
	  }

	 return( TLchgd );

}

const ColorTV::condition ColorTV::coloring(void) {

	Variable::Index ConColor;  // number of consecutive iteration with the same color

	 // green, yellow, red - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 condition lastswing;
	 if( ( -GetdGk() <= -1e-6 ) && ( lastvalue - FiLmbd >=
		       1e-6 * std::max( std::abs( Solver->ReadBestFiVal() ) , 1.0 ) ) )  {
	  lastswing = green;
	  lastgreeniter = NrIter;
	  }
	 else
	  if( ( -GetdGk() > -1e-6 ) && ( lastvalue > FiLmbd ) ) {
	   lastswing = yellow;
	   lastyellowiter = NrIter;
	   }
	  else {
	   lastswing = red;
	   lastrediter = NrIter;
	   }

	 // change beta, if needed - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//	 STPLOG( 1 , std::endl << "           " );

	 switch( lastswing ) {
	  case( green ):  // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   ConColor = NrIter - std::max( lastyellowiter , lastrediter );
//	   STPLOG( 1 , " # green =  " << ConColor << " ~ ScPr = "
//		       << GetdGk() <<  " ~ DeltaFi = " << lastvalue - FiLmbd );
	   if( ( ConColor >= greentestinvl ) && ( Beta < 2 ) ) {
	    lastgreeniter = lastyellowiter = lastrediter = NrIter;
	    Beta *=  2.0;
	    }
	   break;

	  case( yellow ):  //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   ConColor = NrIter - std::max( lastgreeniter , lastrediter );
//	   STPLOG( 2 , " # yellow =  " << ConColor << " ~ ScPr = "
//		       << GetdGk() <<  " ~ DeltaFi = " << lastvalue - FiLmbd );
	   if( ConColor >= yellowtestinvl ) {
	    lastgreeniter = lastyellowiter = lastrediter = NrIter;
	    Beta *= 1.1;
	    }
	   break;

	  case( red ):  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   ConColor = NrIter - std::max( lastgreeniter , lastyellowiter );
//	   STPLOG( 1 , " # red =  " << ConColor << " ~ ScPr = " << GetdGk() );
	   if( ( ConColor >= redtestinvl ) && ( Beta > double( 5e-4 ) ) ) {
	    lastgreeniter = lastyellowiter = lastrediter = NrIter;
	    Beta *= 0.67;
	    }
	   break;

	  }  // end( switch( lastswing ) )- - - - - - - - - - - - - - - - - - - - - -

	 // update the last value - - - - - - - - - - - - - - - - - - - - - - - - - -

	 lastvalue = Solver->ReadFiVal();
	 return( lastswing );



}



/*--------------------------------------------------------------------------*/
/*------------------------- METHODS of FumeroTV ----------------------------*/
/*--------------------------------------------------------------------------*/

FumeroTV::FumeroTV(SubGradSolver* slvr) : Stepsize::Stepsize( slvr  ) {

	this->sigmaMin=stepsizeconf->f_value[par_sigmaMin]; // initial value of beta
    this->r1=stepsizeconf->f_value[par_r1];
	this->BetaZero=stepsizeconf->f_value[par_BetaZero];
	this->eta1=stepsizeconf->f_value[par_eta1];
	this->eta2=stepsizeconf->f_value[par_eta2];




}

void FumeroTV::SetSTPLog(std::ostream* outs, const char lvl) {
	Stepsize::SetSTPLog( outs , lvl );

	 #if( LOG_STP )
	  if( STPLLvl > 1 ) {
	   *STPLog << endl << "FumeroTV: ~ BetaZero = " << BetaZero
		   << " ~ (Eps0 , r1) = (" << sigma << " , " << r1 << ")"
		   << " ~  Eta = {" << eta1 << " , " << eta2 << "}" << endl;
	   }
	 #endif
}

void FumeroTV::Format(void) {
	Stepsize::Format( );
	 Beta = BetaZero;

	 FiLambda = FiRef = Inf<double>();
	 // FiLambda and the reference value are not known
	 LwrBnd = -Inf<double>();  // no lower bound or target level is known

	 etac = rc = 0;
	 r2 = ceil( pow( ( - log( sigmaMin ) / 0.6933 ) , ( 1.0 / 3.26 ) ) * r1 );

	 sigma = 1;
}

void FumeroTV::NewStep(void) {
	// get the function value   - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 FiLambda = ReadFkVal();

	 if( FiLambda == Inf<double>() )
	  std::logic_error( "FumeroTV::GetStepsize: this should not happen" );

	 // try to improve the target level - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 bool TLchgd = UpdateTargetLevel();

	 // update target level- - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( TLchgd )
	  FiLev = sigma * LwrBnd + ( 1 - sigma ) * FiRef;

	 // fix beta - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( Beta < 1e-6 )
	  Beta = 1e-6;

	 // print the level - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 STPLOG( 1 , endl << "           " << " level  = " << -FiLev << " ~  beta = "
		     << Beta << endl << "           " );
}

bool FumeroTV::UpdateTargetLevel(void) {
	bool TLchgd = false;

	 // the target level is updated if any changes of the lower bound occurs - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 Function::FunctionValue LwrBnd_ = 0;//GetOracle()->get_lower_estimate();
	 for(auto a: GetOracle()){
		 LwrBnd_+=a->get_lower_estimate();
	 }
	 if( LwrBnd_ > LwrBnd ) {
	  LwrBnd = LwrBnd_;
	  TLchgd = true;
	  }

	 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( LwrBnd == - Inf<double>() )
	  std::logic_error(
		           "FumeroTV::UpdateTargetLevel: no lower bound is given" );
	 double FiRef_;
	 if( FiRef == Inf<double>() )  // not initialized yet
	  FiRef_ = Inf<double>();
	 else
	  FiRef_ = FiRef - ( 1e-6 * std::max( 1.0 , std::abs( FiRef ) ) );

	 // beta updating- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( FiLambda <= FiRef_ ) {
	  FiRef = FiLambda;
	  if( rc >= r2 )
	   Beta *= 2.0;
	  etac = 0;
	  TLchgd = true;
	  }
	 else  // FiLambda > FiRef_
	  if( rc >= r2 ) {
	   if( ++etac >= eta2 ) {
	    Beta /= 2;
	    etac = 0;
	    }
	   }
	  else
	   if( ++etac >= eta1 ) {
	    Beta = Beta / (2.0 * Beta + 1.0);
	    sigma = ExpFun( ++rc );
	    TLchgd = true;
	    etac = 0;
	    }

	 // failure of the target level updating - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( FiLev < LwrBnd ) {  // change the level
	  FiLev = LwrBnd;
	  TLchgd = true;
	  }

	 return( TLchgd );

}

double FumeroTV::ExpFun(Variable::Index ri) {

	 return( exp( -0.6933 * pow( double(ri) / double( r1 ), 3.26 ) ) );

}




/*--------------------------------------------------------------------------*/
/*------------------------- METHODS of Polyak ------------------------------*/
/*--------------------------------------------------------------------------*/



Polyak::Polyak(SubGradSolver* slvr) : Stepsize::Stepsize( slvr  ) {

	this->Beta=stepsizeconf->f_value[par_Beta];

}

void Polyak::SetSTPLog(std::ostream* outs, const char lvl) {

	Stepsize::SetSTPLog( outs , lvl );

	 #if( LOG_STP )
	  if( STPLLvl > 1 ) {
	   *STPLog << endl << "Polyak: Beta = "<< Beta << " ~ Chi = " <<  LpsFct
		   << endl;
	   }
	 #endif


}

void Polyak::Format(void) {
	MaxBeta = Inf<double>();
	FiLev = - Inf<double>();
}

void Polyak::NewStep(void) {
	// update the target value - - - - - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 UpdateTargetLevel();

	 // some exception  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( FiLev == - Inf<double>() )
	  std::logic_error( "Polyak::GetStepsize: no lower bound is given"  );


}

bool Polyak::UpdateTargetLevel(void) {

	Function::FunctionValue LwrBnd =0;// GetOracle()->get_lower_estimate();
	for(auto a: GetOracle()){
			 LwrBnd+=a->get_lower_estimate();
		 }
	 if( LwrBnd > FiLev ) {
	  FiLev = LwrBnd;
	  return( true );
	  }

	 return( false );

}


