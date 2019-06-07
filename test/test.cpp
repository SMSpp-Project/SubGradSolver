/*--------------------------------------------------------------------------*/
/*-------------------------- File test.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file


 */
/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINES -----------------------------------*/
/*--------------------------------------------------------------------------*/
/*
 */



/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <NDOBlock.h>
#include <SubGradSolver.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>


#include "ThinComputeInterface.h"
#include "Configuration.h"
/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*------------------------------- GLOBALS ----------------------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // reading command line parameters - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


SubGradSolver *solver = new SubGradSolver();


 int pr = solver->get_int_par(1);

 double pr2 =  solver->get_dbl_par(1);

 std::cout << "numero parametri dbl : " << solver->get_num_dbl_par() << std::endl;

 std::cout << "numero parametri int : " <<  solver->get_num_int_par()<<std::endl;

 std::cout << "numero parametri str : " <<  solver->get_num_str_par()<<std::endl;



 for(int i=0;i<solver->get_num_dbl_par();i++){

 std::cout << "numero parametri dbl to str : " <<  solver->dbl_par_idx2str(i) << std::endl;

 }

 for(int i=0;i<solver->get_num_int_par();i++){

  std::cout << "numero parametri int to str : " <<  solver->int_par_idx2str(i) << std::endl;

  }

 //Configuration *d = Configuration::new_Configuration("SimpleConfiguration< std::vector< double > >");

 auto stepconf=Configuration::new_Configuration("SimpleConfiguration< std::vector< double > >");

 auto objstep= dynamic_cast<SimpleConfiguration< std::vector< double > >*>(stepconf);


 std::vector< double > v = {1,1,1,400,10};
 objstep->f_value=v;
 solver->setConfstepsize(stepconf);
 solver->setConfdeflection(stepconf);


// std::cout << "nome : " << solver->name() << std::endl;
//
// //std::ostream myfile("test.txt");
//
//
//
//
// solver->set_log(&std::cout);
//
//
//
//
//
//
//
//bool a=true;
//
//
// std::cout << "compute config name "<<solver->get_ComputeConfig(true)->name() <<std::endl;



// Configuration *d =solver->get_ComputeConfig()->new_Configuration("SimpleConfig_p_p");
// std::cout << "ComputeConfig d  "<<d <<std::endl;
//
// std::cout << "name config anidati "<<d->name() <<std::endl;
//
////
// ComputeConfig *cg1=new ComputeConfig() ;
// std::cout << "ComputeConfig "<<cg1 <<std::endl;
//
// ComputeConfig *cg2=solver->get_ComputeConfig(true,cg1);
// std::cout << "ComputeConfig cg2  "<<cg2 <<std::endl;
//
// //cg2->f_extra_Configuration=d;
//
// std::cout << "ComputeConfig cg2  "<<cg2 <<std::endl;
//
// std::cout << "ComputeConfig d  "<<d <<std::endl;
//
// solver->set_ComputeConfig(cg2);

 //cg1->f_extra_Configuration=d;
// ComputeConfig(cg1);





						// cg->f_extra_Configuration=d;
// cg->f_diff=false;
// solver->set_ComputeConfig(cg);
//
// auto compconf1= solver->get_ComputeConfig(true);
//



 //compconf1->f_extra_Configuration=d;

// solver->set_ComputeConfig(compconf1);

//auto compconf=solver->get_ComputeConfig(true);


//solver->set_ComputeConfig(compconf);




//
// std::cout << solver->get_ComputeConfig(true)->int_pars;
// std::cout << solver->get_ComputeConfig(true)->dbl_pars;
//
// ComputeConfig* solverconf=solver->get_ComputeConfig(true);
//
//
//
// Configuration *c = Configuration::new_Configuration("SimpleConfiguration<int>");
//
// int f=5;
// SimpleConfiguration<int> *sc= new SimpleConfiguration<int>(f);
// SimpleConfiguration<int> *sc1;
// c=sc;
//
// solver->get_ComputeConfig(true)->f_extra_Configuration=c;
//
// //solver->get_ComputeConfig(true)->f_extra_Configuration->serialize("Config_5",eConfigFile);
// //c->serialize("Config_4",eConfigFile);
//
//
// auto d=Configuration::deserialize("Config_4");
//
// sc1=dynamic_cast<SimpleConfiguration<int> *>( d);
//
//
//
//
// std::cout << "f_value " << sc1->f_value ;

//// Configuration *d = Configuration::new_Configuration("SimpleConfiguration<Configuration *,Configuration *>");
//
// ComputeConfig *c1 =new ComputeConfig();
//
// std::pair<std::string,int> v("uno",1);
//
// std::vector<std::pair<std::string, int> > v1(1);//("uno",1);
//
//
// auto cc= c->deserialize("Config_0");
//
//
//
// std::cout <<"c1 conf "<< c1->int_pars << std::endl;
// std::cout <<"c1 name "<< c1->name() << std::endl;
// //std::cout <<"c1 conf "<< c1->int_pars << std::endl;
// auto p= c->new_Configuration("SimpleConfiguration<int>");
//
// solver->get_ComputeConfig(true)->f_extra_Configuration=p;
//
// std::cout  <<"sdfsdfsdfsdf   "<< solver->get_ComputeConfig(true)->dbl_pars;
//
//SimpleConfiguration< std::pair<int,int> > *s= new SimpleConfiguration< std::pair<int,int> >();
//
//s->f_value.first=1;
//
//
//std::filebuf fb;
//fb.open ("test.txt",std::ios::in);
//std::istream os(&fb);
//
//
////p->
// //p << os;
//
// std::cout << "os get" << os.get()<< std::endl;
//
//fb.close();
//
//std::istream *is;
//
//
//
//
//
//std::cout << cc << "c****" << std::endl;
//
////SimpleConfiguration< std::pair<Configuration *,Configuration *> > *sss= new SimpleConfiguration< std::pair<Configuration *,Configuration *> >();
//
////s->f_value.first=1;
//
//
////s->serialize(netCDF::NcFile::ChildrenGrps);
//
//std::cout <<"f_value " << s->f_value << std::endl;
//
//
// std::cout <<"test conf "<< solver->get_ComputeConfig(true,c1)->int_pars << std::endl;
//
// std::cout <<"test conf 2 " << std::endl;
//
//
//solverconf->serialize("Config_0",eConfigFile);
//
//
// std::cout << "nome config  : " <<  p->name() << std::endl;
//
// solverconf->f_extra_Configuration=p;
//
// //p->serialize("Config_0",eConfigFile);
// //auto sp=p->deserialize("Config_0");
//
//
// auto pp=p->new_Configuration("SimpleConfiguration<int>");
//
// std::cout << "f_extra_Configuration   : " << p  << pp<< std::endl;


 Stepsize *stpsz= nullptr;

 Deflection *dfl= nullptr;


//Stepsize *stz=new Stepsize(solver);

 // construct the "original" MCFBlock - - - - - - - - - - - - - - - - - - - -

 //stpsz = dynamic_cast<Stepsize *>( Stepsize( "ColorTV" ) );

 //stpsz= new Stepsize::newStepsize()


//auto s=stpsz->newStepsize("ColorTV",solver);


 //std::cout << "Beta : " << stpsz->GetBeta() << std::endl;


 std::cout << "Beta : " << stpsz->new_Stepsize("ColorTV",solver) << std::endl;

 std::cout << "Beta : " << dfl->new_Deflection("Volume",solver) << std::endl;


//stpsz->Format();
//stpsz->newtest();


  //assert( stpsz );




 //std::cout << stpsz->f_factory();

 //std::cout << stpsz->f_factory();


 //myfile->
 puts("SubGradSolver");

 return( 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
