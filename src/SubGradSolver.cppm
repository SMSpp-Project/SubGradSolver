/*--------------------------------------------------------------------------*/
/*------------------------ File SubGradSolver.cppm -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.subgradsolver, which exports every name
 * that the headers of SubGradSolver declare at namespace scope, together with
 * the modules it depends on; a translation unit can then write
 *
 *     import smspp.subgradsolver;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
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
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "ColorTV.h"
#include "Deflection.h"
#include "FumeroTV.h"
#include "Polyak.h"
#include "PrimalDual.h"
#include "Stepsize.h"
#include "SubGradSolver.h"
#include "Volume.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.subgradsolver;

export import smspp;

export namespace SMSpp_di_unipi_it {

 // ColorTV.h
 using SMSpp_di_unipi_it::ColorTV;

 // Deflection.h
 using SMSpp_di_unipi_it::Deflection;

 // FumeroTV.h
 using SMSpp_di_unipi_it::FumeroTV;

 // Polyak.h
 using SMSpp_di_unipi_it::Polyak;

 // PrimalDual.h
 using SMSpp_di_unipi_it::PrimalDual;

 // Stepsize.h
 using SMSpp_di_unipi_it::Stepsize;

 // SubGradSolver.h
 using SMSpp_di_unipi_it::SubGradSolver;

 // Volume.h
 using SMSpp_di_unipi_it::Volume;
}

/*--------------------------------------------------------------------------*/
/*---------------------- End File SubGradSolver.cppm -----------------------*/
/*--------------------------------------------------------------------------*/
