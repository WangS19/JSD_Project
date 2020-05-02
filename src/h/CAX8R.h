/*****************************************************************************/
/*     To acheive a CAX8R element                                            */
/*     Yaojialun                                                             */
/*     2020/4/14                                                             */
/*****************************************************************************/
#pragma once

#include "Element.h"

using namespace std;

//Plane strain element class

class CAX8R : public CElement
{
private:
	
	//! Calculate the elastic matrix
	double D[4][4];
	void Calculate_D(CAX8RMaterial* material_,double D[][4]);

	double B[4][16];

	//! The number of Gauss point in one direction
	static int NG;

	int NGauss;

	//! Calculate B matrix and Jaccobi matrix
	void Get_BMat(CNode** nodes_, double B[4][16], double& Jac, double& r, double g, double h);
	//For Get_BMat:
	//! Interpolation Function N[8]
	void GetN_2D_SecOrder(double N[8],double g, double h);
	//!dN= [dNdg dNdh]
	void GetdN_2D_SecOrder(double dN[8][2],double g, double h);
	//! Jacobian JMat=[dudg dudh; dwdg dwdh]
	void Get_Jacobian2D(double JMat[2][2], double g, double h, double RZ[2][8]);




public:


	//!	Constructor
	CAX8R();

	//!	Deconstructor
	~CAX8R();

	//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

	//!	Write element data to stream
	virtual void Write(COutputter& output, unsigned int Ele);

	//! Generate location matrix: the global equation number that corresponding to each DOF of the element
	//	Caution:  Equation number is numbered from 1 !
	virtual void GenerateLocationMatrix();

	//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

	//!	Calculate element mass matrix
	virtual void ElementMass(double* stiffness) ;

	//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

	//!	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix(){return 300;};

	//! Get the number of Gauss point in one direction
	virtual void GetNG(int N_G) { NG = N_G; };

	//! Get the number of nodes in one direction
	virtual int GetNN() { return NEN_; };

};
