/*
 * SUBGRADSolver.cpp
 *
 *  Created on: Apr 17, 2019
 *      Author: simone
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/


#include "SUBGRADSolver.h"







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




SUBGRADSolver::SUBGRADSolver() {
	 // initialize algorithmic parameters - - - - - - - - - - - - - - - - - - - -

	KSGPar1=get_int_par(int_par_str2idx("SGPar1"));
	KSGPar2=get_dbl_par(dbl_par_str2idx("SGPar2"));
	KSGPar3=get_int_par(int_par_str2idx("SGPar3"));
	KSGPar4=get_int_par(int_par_str2idx("SGPar4"));
	KSGPar5=get_int_par(int_par_str2idx("SGPar5"));
	srand( KSGPar5 );





	 //if( SGPar2 < 0 )

	//	 this->set_par(this->dbl_par_str2idx("SGPar2"),double(0));


	 // some initializations - - - - - - - - - - - - - - - - - - - - - - - - - - -

//stepsize = nullptr;
	 //deflection = nullptr;
	 //Lambda = LambdaHat = LambdaBest = nullptr;
	 //Gi = dir = nullptr;
	 ub = lb = double (0);

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

SUBGRADSolver::~SUBGRADSolver() {
	// TODO Auto-generated destructor stub
}

void SUBGRADSolver::set_log(std::ostream* log_stream) {
	 Solver::set_log(log_stream);

	 *f_log << std::endl << "SubGrad Solver: Vars = " << "NumVar" << " ("
	 	   << MaxNumVar << ") ~  NItIncr = " << NItIncr <<  " ~ EpsLin = "
	 	   << EpsLin << " ~ FiEps = " << "EInit" << " ~ t* = " << tStar;

	    *f_log << std::endl << " SGParameters = ( " << KSGPar1 << " ~ " << KSGPar2
	 	   << " ~ " << KSGPar3 << " ~ " << ( KSGPar4? "true" : "false")
	 	   <<  " ~ "  << KSGPar5 << " )" << std::endl;


 }

ComputeConfig* SUBGRADSolver::get_ComputeConfig(bool all,
		ComputeConfig* ocfg) const {
	auto tc=ThinComputeInterface::get_ComputeConfig(all,ocfg);

	auto sc = tc->new_Configuration("SimpleConfig_p_p");

	auto s= dynamic_cast<SimpleConfiguration<std::pair< Configuration * , Configuration * > >*>(sc);


	s->f_value.first=this->getConfstepsize();
	s->f_value.second=this->getConfdeflection();


	tc->f_extra_Configuration=s;




	return tc;


}



void SUBGRADSolver::set_ComputeConfig(ComputeConfig* scfg) {

	ThinComputeInterface::set_ComputeConfig(scfg);



}

  Configuration* SUBGRADSolver::getConfdeflection() const {
	return confdeflection;
}

void SUBGRADSolver::setConfdeflection(Configuration*& confdeflection)  {
	this->confdeflection = confdeflection;
}

  Configuration* SUBGRADSolver::getConfstepsize()  const{
	return confstepsize;
}



void SUBGRADSolver::SetStepsize(Stepsize* STP) {

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

//void SUBGRADSolver::SetDeflection(Deflection* Vol) {
//	if( deflection ) {
//	  // if the Stepsize and Deflection objects are actually the same, ensure
//	  // it is not deleted twice
//	  auto StpDef = dynamic_cast<Stepsize *>( deflection );
//	  delete deflection;
//	  if( ! StpDef )
//	   stepsize = nullptr;
//	  }
//
//	 deflection = Vol;
//
//	 if( Oracle )
//	  InitDeflection();
//	 else
//	  alpha = 1;
//}

void SUBGRADSolver::SetFiOracle(C05Function* Fi) {
}

void SUBGRADSolver::SetLambda(std::vector<double> tLambda) {
}

void SUBGRADSolver::KeepBestLambda(const bool KBL) {
}

void SUBGRADSolver::setConfstepsize( Configuration*& confstepsize) {
	this->confstepsize = confstepsize;
}




Function::Vec_FunctionValue SUBGRADSolver::ReadBestSol(const std::vector<Variable::Index> &I , Variable::Index &D ) {

	//I = nullptr;
	D = NumVar;
	return( KpBstL ? LambdaBest : LambdaBar );

}

Function::FunctionValue SUBGRADSolver::ReadBestFiVal(Variable::Index wFi) {

	if( wFi!= Inf<double>() )
	  std::logic_error( "SubGrad::ReadBestFiVal: wFi != INF" );

	 return( FiBest );

}

Function::Vec_FunctionValue SUBGRADSolver::ReadSol(
		const std::vector<Variable::Index> &I , Variable::Index &D) {
	//I = nullptr;
	 D = NumVar;
	 return( LambdaBar );

}

Function::FunctionValue SUBGRADSolver::ReadFiVal( Variable::Index wFi ) {

	if( wFi!= Inf<Variable::Index>() )
	  std::logic_error( "SubGrad::ReadFiVal: wFi != INF"  );

	 return( FiBar );


}

Function::FunctionValue SUBGRADSolver::ReadHatFiVal(void) {

	if( SGPar4 )
	  return( FiHat );
	 else
	  return( Inf<double>() );

}

bool SUBGRADSolver::IsOptimal(Function::FunctionValue eps) {

	double FiL = FiBest;

	 if( FiL == Inf<double>() )
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

}

const Function::FunctionValue* SUBGRADSolver::ReadMult( Variable::Index &I , Variable::Index &D , Variable::Index wFi ) {

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

const  Function::FunctionValue  SUBGRADSolver::ReadLBMult(const Variable::Index wFi) {

	if( wFi!= Inf<Variable::Index>() )
	   std::logic_error( "SubGrad::ReadBestFiVal: wFi != INF" );

	 return( 0 );

}

Variable::Index SUBGRADSolver::NrIter(void) const {
	return( ParIter );

}

const double SUBGRADSolver::ReadFVal( void){

	return FiLambda;

}

void SUBGRADSolver::AllocateUBLB(const Variable::Index strt,
		const Variable::Index num) {
}

double SUBGRADSolver::lbnd(const Variable::Index i) {
	//return( lb ? lb[ i ] : - Inf<double>() );
}

double SUBGRADSolver::ubnd(const Variable::Index i) {
	//return( ub ? ub[ i ] : Inf<double>() );
}

void SUBGRADSolver::EvalGi(Variable::Index wFi) {

	 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( wFi > NrFi )
	  wFi = Inf<Variable::Index>();

	 if( wFi && ! Oracle->get_linearization_constant( wFi ))// NewGi( wFi ) )
	  throw( std::logic_error(
			    "SUBGRADSolver::EvalGi: the subgradient cannot be computed" ) );

	 // fetch the item (subgradient) from the Oracle - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 const std::vector<Variable::Index> SGBse;
	 //Variable::Index SGBDim = Oracle->get_linearization_coefficients() GetGi( Gi , SGBse , wFi? Inf<Variable::Index>() : MaxName );

	 // densify the item  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 //if( SGBse )
	  //Densify( Gi , SGBse , SGBDim , NumVar );

	 // compute the linearization error of Gi[wFi] at Lambda1 - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 Sigma = Oracle->get_value() ;//GetVal();
	 if( KSGPar4 )
	  SigmaHat = Sigma;

	 // write the item  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( wFi && ( ! InnIter ) )
	  Oracle->store_linearization(MaxNCnst);// SetGiName( MaxNCnst );

	 GiEvaltns++;

}

bool SUBGRADSolver::ProjFsb(double Lpt, bool& EmptySet) {
}

void SUBGRADSolver::ProjTng(std::vector<double> Gpt, std::vector<double> Lpt,
		const Variable::Index strt, Variable::Index stp) {
}

void SUBGRADSolver::InitStepsize(void) {
}

void SUBGRADSolver::InitDeflection(void) {
}

void SUBGRADSolver::ChgLambdaHat(void) {
}

void SUBGRADSolver::UpdateSigma(void) {
}

void SUBGRADSolver::UpdtLowerBound(void) {
}

void SUBGRADSolver::MemDealloc(void) {
}
