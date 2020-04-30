/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "LoadCaseData.h"

#include <iomanip>
#include <iostream>

using namespace std;

CLoadCaseData :: ~CLoadCaseData()
{
	delete[] node;
	delete[] dof;
	delete[] load;
	delete[] a_;
	delete[] b_;
	delete[] c_;
	delete[] d_;
	delete[] e_;
}

void CLoadCaseData :: Allocate(unsigned int num)
{
	nloads = num;
	node = new unsigned int[nloads];
	dof = new unsigned int[nloads];
	load = new double[nloads];
	a_ = new double[nloads];
	b_ = new double[nloads];
	c_ = new double[nloads];
	d_ = new double[nloads];
	e_ = new double[nloads];
}; 

//	Read load case data from stream Input
bool CLoadCaseData :: Read(ifstream& Input, unsigned int lcase, int L_type)
{
//	Load case number (LL) and number of concentrated loads in this load case(NL)

	if (L_type == 2)
	{
		return true;
	}
	
	unsigned int LL, NL;

	Input >> LL >> NL;	

	if (LL != lcase + 1) 
	{
		cerr << "*** Error *** Load case must be inputted in order !" << endl 
			 << "   Expected load case : " << lcase + 1 << endl
			 << "   Provided load case : " << LL << endl;

		return false;
	}

	Allocate(NL);

//  L_type	1 -- static load
//			3 -- time-history load
	if (L_type == 1)
	{
		for (unsigned int i = 0; i < NL; i++)
			Input >> node[i] >> dof[i] >> load[i];
	}
	else if (L_type == 3)
	{
		for (unsigned int i = 0; i < NL; i++)
			Input >> node[i] >> dof[i] >> a_[i] >> b_[i] >> c_[i] >> d_[i] >> e_[i];
		Input >> time_knot;
	}
	
	return true;
}

//	Write load case data to stream
void CLoadCaseData::Write(COutputter& output, unsigned int lcase)
{
	for (unsigned int i = 0; i < nloads; i++)
	{
		output << setw(7) << node[i] << setw(13) << dof[i]  << setw(19) << load[i] << endl;
	}
}
