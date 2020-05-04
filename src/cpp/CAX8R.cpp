/*****************************************************************************/
/*     To acheive a plane strain element                                     */
/*     YaoJiaLun                                                             */
/*     2020/4/14                                                             */
/*****************************************************************************/
#pragma once

#include"CAX8R.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#define PI 3.1415926535897932384626433832795
int CAX8R::NG;

//Constructor
CAX8R::CAX8R()
{
	NEN_ = 8;	// Each element has 8 nodes
	nodes_ = new CNode*[NEN_];

	ND_ = 24;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//Deconstructor
CAX8R::~CAX8R()
{
}

void CAX8R::Calculate_D(CAX8RMaterial* material,double D[][4])
{
	double F = material->E / (1.0 + material->mu);
	double G = F / (1.0 - 2.0 * material->mu) * material->mu;
	double H = F + G;
	D[0][0] = H;   D[0][1] = G;   D[0][2] = G;   D[0][3] = 0;
	D[1][0] = G;   D[1][1] = H;   D[1][2] = G;   D[1][3] = 0;
	D[2][0] = G;   D[2][1] = G;   D[2][2] = H;   D[2][3] = 0;
	D[3][0] = 0;   D[3][1] = 0;   D[3][2] = 0;   D[3][3] = F / 2.0;

}

//! Calculate B matrix and Jaccobi matrix
void CAX8R::Get_BMat(CNode** nodes_, double B[4][16], double& Jac, double& r, double g, double h)
{
	//! Construct the shape function
	double N[8];
	GetN_2D_SecOrder(N, g, h);

	//! Calculate the derivative of the shape function
	double dN[8][2];
	GetdN_2D_SecOrder(dN, g, h);

	//! Get NodesCoord (nodes_[k]->XYZ[j])
	double RZ[2][8];
	for (int i = 0; i!=2; i++)
		for (int j = 0; j!=8; j++)
			RZ[i][j] = nodes_[j]->XYZ[i];

	//! Calculate the determinant of Jaccobi matrix
	double JMat[2][2];//JMat=[dudg dwdg; dudh dwdh]
	Get_Jacobian2D(JMat, g, h, RZ);
	
	double det_J ;
	det_J = JMat[0][0] * JMat[1][1] - JMat[1][0] * JMat[0][1];
	Jac = det_J;
	if (det_J < 1e-7) {
		cout << "*** Error *** Some element is singular";
		exit(0);
	}

	//! Calculate the inverse of Jaccobi matrix
	double det_JI = 1.0 / det_J; 
	double JMat_Inv[2][2];//JMat_Inv=[dgdu dhdu; dgdw dhdw]
	JMat_Inv[0][0] = JMat[1][1] * det_JI;
	JMat_Inv[0][1] = -JMat[0][1] * det_JI;
	JMat_Inv[1][0] = -JMat[1][0] * det_JI;
	JMat_Inv[1][1] = JMat[0][0] * det_JI;

	//! Calculate the B matrix
	r = 0;
	for (int i = 0; i!=8; i++)
		r += N[i]*RZ[0][i];
	
	for(int i = 0; i!=8; i++)
	{
		// B[0][8] = [dN1dr 0 ... dN8dr 0]
		// dNdu = dNdg*dgdu+dNdh*dhdu = JMat_Inv[0][0]*dN[i][0]+JMat_Inv[0][1]*dN[i][1]
		B[0][i*2]   = JMat_Inv[0][0]*dN[i][0]+JMat_Inv[0][1]*dN[i][1];
		B[0][i*2+1] = 0;

		// B[1][8] = [0 dN1dz ... 0 dN8dz]
		// dNdw = dNdg*dgdw+dNdh*dhdw = JMat_Inv[1][0]*dN[i][0]+JMat_Inv[1][1]*dN[i][1]
		B[1][i*2]   =  0;
		B[1][i*2+1] =  JMat_Inv[1][0]*dN[i][0]+JMat_Inv[1][1]*dN[i][1];

		// B[2][8] = [N1/r 0 ... N8/r 0]
		B[2][i*2]   = N[i]/r;
		B[2][i*2+1] = 0;

		// B[3][8] = [dN1dw dN1du ... dN8dw dN8du]
		B[3][i*2]   =  B[1][i*2+1];
		B[3][i*2+1] =  B[0][i*2];
	}


}
void CAX8R::GetN_2D_SecOrder(double N[8],double g, double h)
{
	double gp1 = 1+g;
	double gm1 = 1-g;
	double hp1 = 1+h;
	double hm1 = 1-h;
	double gp1ph = 1+g+h;
	double gp1mh = 1+g-h;
	double gm1ph = 1-g+h;
	double gm1mh = 1-g-h;
	
	N[0] =  -0.25*gm1*hm1*gp1ph;
	N[1] =  -0.25*gp1*hm1*gm1ph;
	N[2] =  -0.25*gp1*hp1*gm1mh;
	N[3] =  -0.25*gm1*hp1*gp1mh;
	N[4] =  0.5*gm1*gp1*hm1;
	N[5] =  0.5*hm1*hp1*gp1;
	N[6] =  0.5*gm1*gp1*hp1;
	N[7] =  0.5*hm1*hp1*gm1;
}
void CAX8R::GetdN_2D_SecOrder(double dN[8][2],double g, double h)
{
	double gp1 = 1+g;
	double gm1 = 1-g;
	double hp1 = 1+h;
	double hm1 = 1-h;
	double gp1ph = 1+g+h;
	double gp1mh = 1+g-h;
	double gm1ph = 1-g+h;
	double gm1mh = 1-g-h;
	//dNdg
	dN[0][0] =  0.25*hm1*gp1ph-0.25*gm1*hm1;			//-0.25*gm1*hm1*gp1ph;
	dN[1][0] =  -0.25*hm1*gm1ph+0.25*gp1*hm1;			//-0.25*gp1*hm1*gm1ph;
	dN[2][0] =  -0.25*hp1*gm1mh+0.25*gp1*hp1;			//-0.25*gp1*hp1*gm1mh;
	dN[3][0] =  0.25*hp1*gp1mh-0.25*gm1*hp1;			//-0.25*gm1*hp1*gp1mh;
	dN[4][0] =  -0.5*gp1*hm1+0.5*gm1*hm1;				//0.5*gm1*gp1*hm1;
	dN[5][0] =  0.5*hm1*hp1;							//0.5*hm1*hp1*gp1;
	dN[6][0] =  -0.5*gp1*hp1+0.5*gm1*hp1;				//0.5*gm1*gp1*hp1;
	dN[7][0] =  -0.5*hm1*hp1;							//0.5*hm1*hp1*gm1;
	//dNdh
	dN[0][1] =  0.25*gm1*gp1ph-0.25*gm1*hm1;			//-0.25*gm1*hm1*gp1ph;
	dN[1][1] =  0.25*gp1*gm1ph-0.25*gp1*hm1;			//-0.25*gp1*hm1*gm1ph;
	dN[2][1] =  -0.25*gp1*gm1mh+0.25*gp1*hp1;			//-0.25*gp1*hp1*gm1mh;
	dN[3][1] =  -0.25*gm1*gp1mh+0.25*gm1*hp1;			//-0.25*gm1*hp1*gp1mh;
	dN[4][1] =  -0.5*gm1*gp1;							//0.5*gm1*gp1*hm1;
	dN[5][1] =  -0.5*hp1*gp1+0.5*hm1*gp1;				//0.5*hm1*hp1*gp1;
	dN[6][1] =  0.5*gm1*gp1;							//0.5*gm1*gp1*hp1;
	dN[7][1] =  -0.5*hp1*gm1+0.5*hm1*gm1;				//0.5*hm1*hp1*gm1;
}
void CAX8R::Get_Jacobian2D(double JMat[2][2], double g, double h, double RZ[2][8])
{
	//JMat=[dudg dwdg; dudh dwdh]  RZ= [u[8]; w[8]]
	double dN[8][2];
	GetdN_2D_SecOrder(dN, g, h);
	for (int i = 0; i!=2; i++)
	{
		for (int j = 0; j!=2; j++)
		{
			JMat[i][j] = 0;
			for (int k = 0; k!=8; k++)
				JMat[i][j] += dN[k][i]*RZ[j][k];
		}
	}
	
}
//public
void CAX8R::ElementStiffness(double* Matrix)
{
	int n = SizeOfStiffnessMatrix();
	clear(Matrix, n);

	CAX8RMaterial* material_ = dynamic_cast<CAX8RMaterial*>(ElementMaterial_);
	//Calculate the elastric matrix
	Calculate_D(material_,D);

	double KK[16][16];
	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++)
			KK[i][j] = 0;
	}
	
	double SS[4]; //! the column of SS = D * B
	double r;
	for (int Gx = 0; Gx < NG; Gx++) {
		double XI = XG[Gx][NG - 1];
		for (int Gy = 0; Gy < NG; Gy++) {
			double YI = XG[Gy][NG - 1];
			double WT = WGT[Gx][NG-1] * WGT[Gy][NG-1];
			//! Calculate B matrix and Jaccobi matrix
			Get_BMat(nodes_, B, det_J, r, XI, YI);

			//! Calculate the element stiffness matrix
			for (int j = 0; j < 16; j++) {
				for (int k = 0; k < 4; k++) {
					SS[k] = 0.0;
					for (int l = 0; l < 4; l++) {
						SS[k] += D[k][l] * B[l][j];
					}
				}
				for (int i = j; i < 16; i++) {
					double STIFF = 0.0;
					for (int l = 0; l < 4; l++) {
						STIFF += B[l][i] * SS[l];
					}
					KK[j][i] += 2*PI*r*STIFF * WT * det_J;
				}
			}
		}
	}

	//! Extend KK into 24*24 matrix
	double K[24][24];
	for (int i = 0; i < 24; i++) {
		for (int j = 0; j < 24; j++) {
			K[i][j] = 0.0;
		}
	}
	int m1, m2;
	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) {
			int k1 = i / 2;
			if (i % 2 == 0)
				m1 = 0;
			else
				m1 = 1;
			int k2 = j / 2;
			if (j % 2 == 0)
				m2 = 0;
			else
				m2 = 1;
			K[3 * k1 + m1][3 * k2 + m2] = KK[i][j];

		}
	}

	//! Store K into Matrix
	int k = 0;
	for (int i = 0; i < 24; i++) {
		for (int j = i; j >= 0; j--) {
			Matrix[k] = K[j][i];
			k += 1;
		}
	}
}
void CAX8R::ElementStress(double* stress/* 4*4 */, double* Displacement)
{
	CAX8RMaterial* material_ = dynamic_cast<CAX8RMaterial*>(ElementMaterial_);
	Calculate_D(material_, D);
	// Make up the whole displacement vector
	double dis[16];
	double r;
	int k = 0;
	for (int i = 0; i < NEN_; i++) {
		for (unsigned int j = 0; j < 2; j++)
		{
			if (nodes_[i]->bcode[j] == 0)
			{
				dis[k] = 0;
			}
			else
			{
				dis[k] = Displacement[nodes_[i]->bcode[j] - 1];
			}
			k += 1;
		}
	}
	k = 0;
	double stan[4], strs[4]; // The strain and stress
	for (int Gx = 0; Gx < NG; Gx++) {
		double XI = XG[Gx][NG - 1];
		for (int Gy = 0; Gy < NG; Gy++) {
			for (int i = 0; i < 4; i++) {
				stan[i] = 0.0;
				strs[i] = 0.0;
			}
			double YI = XG[Gy][NG - 1];
			//! Calculate B matrix and Jaccobi matrix
			Get_BMat(nodes_, B, det_J, r, XI, YI);
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 16; j++) {
					stan[i] += B[i][j] * dis[j];
				}
			}
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					strs[i] += D[i][j] * stan[j];
				}
			}
			for (int i = 0; i < 4; i++) {
				stress[4 * k + i] = strs[i];
			}
			k += 1;
		}
	}
}
bool CAX8R::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N;	// element number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl
			<< "    Expected element : " << Ele + 1 << endl
			<< "    Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4, N5, N6, N7, N8;	// The node number in anticlockwise

	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8>> MSet;
	ElementMaterial_ = dynamic_cast<CAX8RMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	nodes_[4] = &NodeList[N5 - 1];
	nodes_[5] = &NodeList[N6 - 1];
	nodes_[6] = &NodeList[N7 - 1];
	nodes_[7] = &NodeList[N8 - 1];

	return true;
}

void CAX8R::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) 
		<< nodes_[0]->NodeNumber << setw(9) << nodes_[1]->NodeNumber << setw(9) 
		<< nodes_[2]->NodeNumber << setw(9) << nodes_[3]->NodeNumber << setw(9) 
		<< nodes_[4]->NodeNumber << setw(9) << nodes_[5]->NodeNumber << setw(9) 
		<< nodes_[6]->NodeNumber << setw(9) << nodes_[7]->NodeNumber << setw(9) 
		<< setw(12) << ElementMaterial_->nset << endl;
}

void CAX8R::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 0; D < 3; D++)
			LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

void CAX8R::ElementMass(double* Matrix) 
{
	int n = SizeOfStiffnessMatrix();
	clear(Matrix, n);

	CAX8RMaterial* material_ = dynamic_cast<CAX8RMaterial*>(ElementMaterial_);
	double ROU = material_->rou;

	double MM[16][16];
	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++)
			MM[i][j] = 0;
	}
	
	
	double r;
	for (int Gx = 0; Gx < NG+1; Gx++) //质量阵采用3*3积分才与ABAQUS一致所以NG+1
	{
		double g = XG[Gx][NG];       //质量阵采用3*3积分才与ABAQUS一致所以后面的指标为NG，不是NG-1
		for (int Gy = 0; Gy < NG+1; Gy++) 
		{
			double h = XG[Gy][NG];
			double WT = WGT[Gx][NG] * WGT[Gy][NG];

			//! Construct the shape function
			double N[8];
			GetN_2D_SecOrder(N, g, h);

			//! Get NodesCoord (nodes_[k]->XYZ[j])
			double RZ[2][8];
			for (int i = 0; i!=2; i++)
				for (int j = 0; j!=8; j++)
					RZ[i][j] = nodes_[j]->XYZ[i];

			// Cal r
			r = 0;
			for (int i = 0; i!=8; i++)
				r += N[i]*RZ[0][i];

			//! Calculate the determinant of Jaccobi matrix
			double JMat[2][2];//JMat=[dudg dwdg; dudh dwdh]
			Get_Jacobian2D(JMat, g, h, RZ);

			double det_Jac ;
			det_Jac = JMat[0][0] * JMat[1][1] - JMat[1][0] * JMat[0][1];
			if (det_Jac < 1e-7) {
				cout << "*** Error *** Some element is singular";
				exit(0);
			}

			double scale = ROU*2*PI*r*det_Jac*WT;
			//! Calculate MM=scale*rou*NT*N    N[2*16] = [N1 0...N8 0; 0 N1...0 N8]
			//  NTN[16*16] = [N1N1 0 N1N2 0 ... N1N8 0;
            //               0 N1N1 0 N1N2 ... 0 N1N8;
			//               ...;
			//               N8N1 0 ... N8N8 0;
			//               0 N8N1 ... 0 N8N8; ]
			for (unsigned int m = 0; m!=8; m++)
			{
				for (unsigned int n = 0; n!=8; n++)
				{
					MM[m*2][n*2] += scale*N[m]*N[n];
				}
					
			}
		}
	}
	for (unsigned int m = 0; m!=8; m++)
	{
		for (unsigned int n = 0; n!=8; n++)
		{
			MM[m*2+1][n*2+1] = MM[m*2][n*2];
		}
					
	}
	// CAX8R单元质量阵用
	//对角元放大，系数1.1842满足质量守恒,转化为集中质量阵
	for (int i = 0; i < 16; i++) {
		for(int j = 0; j < 16; j++){
			if(i==j)
				MM[i][j] = 1.1842*MM[i][j];
			else
				MM[i][j] = 0;
		}
		
	}
	//! Extend MM into 24*24 matrix M
	double M[24][24];
	for (int i = 0; i < 24; i++) {
		for (int j = 0; j < 24; j++) {
			M[i][j] = 0.0;
		}
	}
	int m1, m2;
	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) {
			int k1 = i / 2;
			if (i % 2 == 0)
				m1 = 0;
			else
				m1 = 1;
			int k2 = j / 2;
			if (j % 2 == 0)
				m2 = 0;
			else
				m2 = 1;
			M[3 * k1 + m1][3 * k2 + m2] = MM[i][j];

		}
	}

	//! Store M into Matrix
	int k = 0;
	for (int i = 0; i < 24; i++) {
		for (int j = i; j >= 0; j--) {
			Matrix[k] = M[j][i];
			k += 1;
		}
	}
}