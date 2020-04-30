/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"

using namespace std;

//! Class LoadData is used to store load data
class CLoadCaseData
{
public:

	unsigned int nloads;	//!< Number of concentrated loads in this load case
	unsigned int* node;		//!< Node number to which this load is applied
	unsigned int* dof;		//!< Degree of freedom number for this load component
	double* load;			//!< Magnitude of load
	double time_knot;		//!< Store the time knot of the load
	double* a_;				//!< Store the function of time history load Q = a_ + b_*t + c_ * sin(d_*t + e_)
	double* b_;
	double* c_;
	double* d_;
	double* e_;

public:

	CLoadCaseData() : nloads(0), node(NULL), dof(NULL), load(NULL), time_knot(NULL), a_(NULL), b_(NULL), c_(NULL), d_(NULL), e_(NULL){};
	~CLoadCaseData();

//!	Set nloads, and new array node, dof and load
	void Allocate(unsigned int num);

//!	Read load case data from stream Input
//  Load_type   1 -- static load   2 -- no read   3 -- load varies with time
	bool Read(ifstream& Input, unsigned int lcase, int Load_type);

//!	Write load case data to stream
	void Write(COutputter& output, unsigned int lcase);
};
