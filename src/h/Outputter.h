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

#include <fstream>
#include <iostream>

using namespace std;

//! Outputer class is used to output results
class COutputter
{
private:

//!	File stream for output
	ofstream OutputFile;

//! Enlarge the dis by a scale
	double Dis_scale ;//= 100.0;

protected:

//!	Constructor
	COutputter(string FileName);

//!	Designed as a single instance class
	static COutputter* _instance;
//! Another instance class for tecplot
	static COutputter* tec_instance;
//! Another instance class for history message
	static COutputter* his_instance;

public:

//!	Return pointer to the output file stream
	inline ofstream* GetOutputFile() { return &OutputFile; }

//!	Return the single instance of the class
	static COutputter* Instance(string FileName = " ");
//! Return the tecplot instance of the class
	static COutputter* Tec_Instance(string FileName = " ");
//! Return the history instance of the class
	static COutputter* His_Instance(string FileName = " ");

//!	Output current time and date
	void PrintTime(const struct tm * ptm, COutputter& output);

//!	Output logo and heading 
	void OutputHeading();

//!	Output nodal point data
	void OutputNodeInfo();

//!	Output equation numbers
	void OutputEquationNumber();

//!	Output element data
	void OutputElementInfo();

//!	Output bar element data
	void PrintBarElementData(unsigned int EleGrp);

//!	Output Q4 element data
	void PrintQ4ElementData(unsigned int EleGrp);

//!	Output CAX8R element data
	void PrintAX8RElementData(unsigned int EleGrp);//yjl

//!	Output load data 
	void OutputLoadInfo(); 

//!	Output displacement data
	void OutputNodalDisplacement(unsigned int lcase);

//!	Output element stresses 
	void OutputElementStress();

//!	Print total system data
	void OutputTotalSystemData();

//! Output into tecplot (for stastic and modal analysis)
	void OutputTecplot(int step);

//! Overload: Output into tecplot (for dynamics analysis)
	void OutputTecplot(double time, double* dis);

//! Print motion of certain freedom against time
	void OutputHisMessage(double time, double* dis, double* vel, double* acc, int N, int* Freedoms);

//! Overload the operator <<
	template <typename T>
	COutputter& operator<<(const T& item) 
	{
		std::cout << item;
		OutputFile << item;
		return *this;
	}

	typedef std::basic_ostream<char, std::char_traits<char> > CharOstream;
	COutputter& operator<<(CharOstream& (*op)(CharOstream&)) 
	{
		op(std::cout);
		op(OutputFile);
		return *this;
	}

#ifdef _DEBUG_

//!	Print banded and full stiffness matrix for debuging
	void PrintStiffnessMatrix();

//!	Print banded and full mass matrix for debuging
	void PrintMassMatrix();

//!	Print address of diagonal elements for debuging
	void PrintDiagonalAddress();

//!	Print column heights for debuging
	void PrintColumnHeights();

//!	Print displacement vector for debuging
	void PrintDisplacement(unsigned int loadcase);

#endif

};
