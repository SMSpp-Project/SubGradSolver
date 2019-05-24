/*
 * LuckFi.h
 *
 *  Created on: May 22, 2019
 *      Author: simone
 */


/*--------------------------------------------------------------------------*/
/*---------------------------- File LuckFi.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Definition of the LuckFi class, which implements the C05Function interface
 * and provides a set of 29 test problems for nonsmmoth unconstrained
 * optimization. The description of test problems is reported by
 * Ladislav Luksan and Jan Vlcek, it is available at
 *
 *     http://www.cs.cas.cz/~luksan/test.html
 *
 * and in
 *
 *     "Test problems for nonsmooth unconstrained optimization and
 *      linearity constrained optimization"
 *     Technical Report 798, Institute of Computer Science,
 *     Academy of Science of the Czech Republic
 *
 * \version 1.10
 *
 * \date 27 - 01 - 2019
 *
 * \author  Annabella Astorino \n
 *          Istituto di Calcolo e Reti ad Alte Prestazioni - CNR \n
 *
 * \author  Antonio Frangioni \n
 *          Operations Research Group \n
 *          Dipartimento di Informatica \n
 *          Università di Pisa \n
 *
 * \author  Manlio Gaudioso \n
 *          Dipartimento di Elettronica Informatica e Sistemistica \n
 *          Universita' della Calabria \n
 *
 * \author  Enrico Gorgone \n
 *          Dipartimento di Elettronica Informatica e Sistemistica \n
 *          Universita' della Calabria \n
 *
 * Copyright &copy 2009 - 2019 by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef LUCKFI_H_
#define LUCKFI_H_


/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/


#include "C05Function.h"



/*--------------------------------------------------------------------------*/
/*------------------------ NAMESPACE and USINGS ----------------------------*/
/*--------------------------------------------------------------------------*/


using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS LuckFi ------------------------------*/
/*--------------------------------------------------------------------------*/


class LuckFi :public C05Function  {

	/*--------------------------------------------------------------------------*/
	/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
	/*--------------------------------------------------------------------------*/



public:


	enum int_par_type_LF{

		NameF=intLastParC0F, //[1]: the function that the FiOracle returns, chosen in the list
        					//of available 25 ones below

		NrVr, // [2]: a few functions of the list, namely 26, 27 and 28, have a
				// parametric number of variable, that is fixed to NrVr

		NrCmp, 	// [2]: function 28 (MaxQR) has a parametric number of components,
        		//i.e., it is the sum of NrCmp functions

		seed,  // [1]: function 28 (MaxQR) is randomly generated, with this seed

	};


	/*--------------------------------------------------------------------------*/
	/*--------------------------- PUBLIC METHODS -------------------------------*/
	/*--------------------------------------------------------------------------*/
	/*---------------------------- CONSTRUCTOR ---------------------------------*/
	/*--------------------------------------------------------------------------*/




	LuckFi();
	virtual ~LuckFi();
};

#endif /* LUCKFI_H_ */
