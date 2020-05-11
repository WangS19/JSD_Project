/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> Area >> rho;	// Young's modulus and section area and density

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << Area << endl;
}

//	Read material data from stream Input
bool CQ4Material::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl
			<< "    Expected set : " << mset + 1 << endl
			<< "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> miu >> rho >> t >> ss;	// Young's modulus and Poisson ratio and thickness and ss(0-plane strain;1-plane stress)

	return true;
}

//	Write material data to Stream
void CQ4Material::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset + 1 << setw(16) << E << setw(16) << miu << setw(16) << t << endl;
}

// CAX8R yjl
bool CAX8RMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl
			<< "    Expected set : " << mset + 1 << endl
			<< "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> mu >> rou;	

	return true;
}
bool CAX8RMaterial::ReadInp(ifstream& Input, unsigned int mset, streampos pMaterial)
{
	nset = mset + 1 ;
	Input.seekg( pMaterial, ios_base :: beg);
	string line;
	char *c = new char [8];
	memset(c,0,9);
	Input.read(c,8);
	while(!Input.eof())
	{
		Input.seekg( -8, ios :: cur);
		getline(Input,line);
		size_t comma = 0;
		if(!(stricmp(c,"*Density")))
		{
			getline(Input,line);
			comma = line.find (',',0);
			rou = atof(line.substr(0,comma).c_str());
		}
		else if(!(stricmp(c,"*Elastic")))
		{
			getline(Input,line);
			comma = line.find (',',0);
			E = atof(line.substr(0,comma).c_str());
			mu = atof(line.substr(comma + 1 , line.length() - comma - 1).c_str());
			break;
		}
		Input.read(c,8);
	}
	if( Input.eof())
	{
		Input.clear();  //clear eof!!!
		Input.seekg(0, ios_base::beg);
		return false;
	}

	return true;
}

bool CQ4Material :: ReadInp(ifstream& Input, unsigned int mset, streampos pMaterial)
{
	nset = mset + 1 ;
	Input.seekg( pMaterial, ios_base :: beg);
	string line;
	char *c = new char [8];
	memset(c,0,9);
	Input.read(c,8);
	while(!Input.eof())
	{
		Input.seekg( -8, ios :: cur);
		getline(Input,line);
		size_t comma = 0;
		if(!(stricmp(c,"*Density")))
		{
			getline(Input,line);
			comma = line.find (',',0);
			rho = atof(line.substr(0,comma).c_str());
		}
		else if(!(stricmp(c,"*Elastic")))
		{
			getline(Input,line);
			comma = line.find (',',0);
			E = atof(line.substr(0,comma).c_str());
			miu = atof(line.substr(comma + 1 , line.length() - comma - 1).c_str());
			break;
		}
		Input.read(c,8);
	}
	if( Input.eof())
	{
		Input.clear();  //clear eof!!!
		Input.seekg(0, ios_base::beg);
		return false;
	}

	//t & ss to be clear !!!
	Input.seekg(0, ios_base::beg);
	Input.read(c,8);
	while(Input.tellg() < pMaterial)
	{
		Input.seekg(-8, ios :: cur);
		getline(Input,line);
		if (!stricmp(c,"*Solid S"))
			break;
		Input.read(c,8);
	}
	size_t comma = 0;
	getline(Input,line);
	comma = line.find(',',0);
	t = atof(line.substr(0,comma).c_str());
	ss = 1;

	return true;
}
bool CBarMaterial :: ReadInp(ifstream& Input, unsigned int mset, streampos pMaterial)
{
	nset = mset + 1 ;
	Input.seekg( pMaterial, ios_base :: beg);
	string line;
	char *c = new char [8];
	memset(c,0,9);
	Input.read(c,8);
	while(!Input.eof())
	{
		Input.seekg( -8, ios :: cur);
		getline(Input,line);
		size_t comma = 0;
		if(!(stricmp(c,"*Density")))
		{
			getline(Input,line);
			comma = line.find (',',0);
			rho = atof(line.substr(0,comma).c_str());
		}
		else if(!(stricmp(c,"*Elastic")))
		{
			getline(Input,line);
			comma = line.find (',',0);
			E = atof(line.substr(0,comma).c_str());
	//		mu = atof(line.substr(comma + 1 , line.length() - comma - 1).c_str());
			break;
		}
		Input.read(c,8);
	}
	if( Input.eof())
	{
		Input.clear();  //clear eof!!!
		Input.seekg(0, ios_base::beg);
		return false;
	}
	Input.seekg(0, ios_base::beg);
	Input.read(c,8);
	while(Input.tellg() < pMaterial)
	{
		Input.seekg(-8, ios :: cur);
		getline(Input,line);
		if (!stricmp(c,"*Solid S"))
			break;
		Input.read(c,8);
	}
	size_t comma = 0;
	getline(Input,line);
	comma = line.find(',',0);
	Area = atof(line.substr(0,comma).c_str());
	return true;
}
void CAX8RMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset + 1 << setw(16) << E << setw(16) << mu << setw(16) << rou << endl;
}