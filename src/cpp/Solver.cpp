/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Solver.h"

#include <cmath>
#include <cfloat>
#include <limits.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

//	Clear an array
template <class type> void clear(type* a, unsigned int N)
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CSolver::CSolver(CSkylineMatrix<double>* K) : K(K) {};

CSolver::CSolver(CSkylineMatrix<double>* K, CSkylineMatrix<double>* M) : K(K), M(M) {};	

// Calculate the lumped mass matrix
void CSolver::Cal_LM()
{
	unsigned int* DiagonalAddress = M->GetDiagonalAddress();
	int NEQ = M->dim();
	L_M = new double[NEQ];
	clear(L_M, NEQ);
	for (int i = 0; i < NEQ; i++) {
		for (int j = 0; j < NEQ; j++) {
			int H;
			if (j >= i)
				H = DiagonalAddress[j + 1] - DiagonalAddress[j];
			else
				H = DiagonalAddress[i + 1] - DiagonalAddress[i];
			if (j - i - H >= 0 || i - j - H >= 0)
			{
				L_M[i] += 0.0;
			}
			else
			{
				L_M[i] += (*M)(i + 1, j + 1);
			}

		}
	}
}


// Perform L*D*L(T) factorization of any of the matrix with skyline form
void CSolver::LDLT(CSkylineMatrix<double>* KK)
{
	unsigned int N = KK->dim();
	unsigned int* ColumnHeights = KK->GetColumnHeights();   // Column Hights

	for (unsigned int j = 2; j <= N; j++)      // Loop for column 2:n (Numbering starting from 1)
	{
		// Row number of the first non-zero element in column j (Numbering starting from 1)
		unsigned int mj = j - ColumnHeights[j - 1];

		for (unsigned int i = mj + 1; i <= j - 1; i++)	// Loop for mj+1:j-1 (Numbering starting from 1)
		{
			// Row number of the first nonzero element in column i (Numbering starting from 1)
			unsigned int mi = i - ColumnHeights[i - 1];

			double C = 0.0;
			for (unsigned int r = max(mi, mj); r <= i - 1; r++)
				C += (*KK)(r, i) * (*KK)(r, j);		// C += L_ri * U_rj

			(*KK)(i, j) -= C;	// U_ij = K_ij - C
		}

		for (unsigned int r = mj; r <= j - 1; r++)	// Loop for mj:j-1 (column j)
		{
			double Lrj = (*KK)(r, j) / (*KK)(r, r);	// L_rj = U_rj / D_rr
			(*KK)(j, j) -= Lrj * (*KK)(r, j);	// D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
			(*KK)(r, j) = Lrj;
		}

		if (fabs((*KK)(j, j)) <= FLT_MIN)
		{
			cerr << "*** Error *** Stiffness matrix is not positive definite !" << endl
				<< "    Euqation no = " << j << endl
				<< "    Pivot = " << (*KK)(j, j) << endl;

			exit(4);
		}
	}
}

//! Back substitute
void CSolver::BackSubstitution(CSkylineMatrix<double>* KK, double* Force)
{
	{
		unsigned int N = KK->dim();
		unsigned int* ColumnHeights = KK->GetColumnHeights();   // Column Hights

															   //	Reduce right-hand-side load vector (LV = R)
		for (unsigned int i = 2; i <= N; i++)	// Loop for i=2:N (Numering starting from 1)
		{
			unsigned int mi = i - ColumnHeights[i - 1];

			for (unsigned int j = mi; j <= i - 1; j++)	// Loop for j=mi:i-1
				Force[i - 1] -= (*KK)(j, i) * Force[j - 1];	// V_i = R_i - sum_j (L_ji V_j)
		}

		//	Back substitute (Vbar = D^(-1) V, L^T a = Vbar)
		for (unsigned int i = 1; i <= N; i++)	// Loop for i=1:N
			Force[i - 1] /= (*KK)(i, i);	// Vbar = D^(-1) V

		for (unsigned int j = N; j >= 2; j--)	// Loop for j=N:2
		{
			unsigned int mj = j - ColumnHeights[j - 1];

			for (unsigned int i = mj; i <= j - 1; i++)	// Loop for i=mj:j-1
				Force[i - 1] -= (*KK)(i, j) * Force[j - 1];	// a_i = Vbar_i - sum_j(L_ij Vbar_j)
		}
	};
}

// LDLT facterization
void CLDLTSolver::LDLT()
{
	unsigned int N = K->dim();
    unsigned int* ColumnHeights = K->GetColumnHeights();   // Column Hights

	for (unsigned int j = 2; j <= N; j++)      // Loop for column 2:n (Numbering starting from 1)
	{
        // Row number of the first non-zero element in column j (Numbering starting from 1)
		unsigned int mj = j - ColumnHeights[j-1];
        
		for (unsigned int i = mj+1; i <= j-1; i++)	// Loop for mj+1:j-1 (Numbering starting from 1)
		{
            // Row number of the first nonzero element in column i (Numbering starting from 1)
			unsigned int mi = i - ColumnHeights[i-1];

			double C = 0.0;
			for (unsigned int r = max(mi, mj); r <= i-1; r++)
				C += (*K)(r,i) * (*K)(r,j);		// C += L_ri * U_rj

			(*K)(i,j) -= C;	// U_ij = K_ij - C
		}

		for (unsigned int r = mj; r <= j-1; r++)	// Loop for mj:j-1 (column j)
		{
			double Lrj = (*K)(r,j) / (*K)(r,r);	// L_rj = U_rj / D_rr
			(*K)(j,j) -= Lrj * (*K)(r,j);	// D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
			(*K)(r,j) = Lrj;
		}

        if (fabs((*K)(j,j)) <= FLT_MIN)
        {
            cerr << "*** Error *** Stiffness matrix is not positive definite !" << endl
            	 << "    Euqation no = " << j << endl
            	 << "    Pivot = " << (*K)(j,j) << endl;
            
            exit(4);
        }
    }
};

// Solve displacement by back substitution
void CLDLTSolver::BackSubstitution(double* Force)
{
	unsigned int N = K->dim();
    unsigned int* ColumnHeights = K->GetColumnHeights();   // Column Hights

//	Reduce right-hand-side load vector (LV = R)
	for (unsigned int i = 2; i <= N; i++)	// Loop for i=2:N (Numering starting from 1)
	{
        unsigned int mi = i - ColumnHeights[i-1];

		for (unsigned int j = mi; j <= i-1; j++)	// Loop for j=mi:i-1
			Force[i-1] -= (*K)(j,i) * Force[j-1];	// V_i = R_i - sum_j (L_ji V_j)
	}

//	Back substitute (Vbar = D^(-1) V, L^T a = Vbar)
	for (unsigned int i = 1; i <= N; i++)	// Loop for i=1:N
		Force[i-1] /= (*K)(i,i);	// Vbar = D^(-1) V

	for (unsigned int j = N; j >= 2; j--)	// Loop for j=N:2
	{
        unsigned int mj = j - ColumnHeights[j-1];

		for (unsigned int i = mj; i <= j-1; i++)	// Loop for i=mj:j-1
			Force[i-1] -= (*K)(i,j) * Force[j-1];	// a_i = Vbar_i - sum_j(L_ij Vbar_j)
	}
};



//! To solve the generalized eigenproblem with the Lanczos method
void CModal::Lanczos()					
{
	int NEQ = K->dim();
	unsigned int* MAXA = K->GetDiagonalAddress();
	int NWK = K->size();
	double* KK = new double[NWK];
	for (int i = 0; i < NWK; i++){
		KK[i] = (*K)(i+1);
	}
	double M[5] = { 1.0, 1.0, 1.0, 1.0, 0.5 };
	int NWM = 5;

	Eig_value = new double[Nroot];
	Eig_vector = new vector<double>[Nroot];
	for (int i = 0; i < Nroot; i++) {
		Eig_value[i] = 0.0;
		Eig_vector[i].resize(NEQ);
	}

	N_iv = min(Nroot * 2, Nroot + 8);
	N_iv = min(N_iv, NEQ);
	int N_p1 = N_iv + 1;
	int N_conv = 0;		//Number of eigen vectors converged
	int N_ite = 0;		//Number of iteration time
	double shift = 0.0;	//The shift value

	//! === Loop for restart ===
	while (N_conv < Nroot && N_ite <= N_ite_max)
	{
		N_ite += 1;

		//! === Set the shift ===
		if (N_ite > 1 && N_conv > 0) {
			double shift = Eig_value[N_conv - 1];

			if (NWM == NEQ) {	// The lumped mass matrix is used
				for (int i = 0; i < NWM; i++) {
					int j = MAXA[i];
					KK[j] -= shift * M[i];
				}
			}
			else {		// The consistent mass matrix is used
				for (int i = 0; i < NWM; i++) {
					KK[i] -= shift * M[i];
				}
			}
		}

	}

}

// To solve the standard eigenvalue problem with the Jacobi method
void CModal::Jacobi()
{

}

// Use the Gram-Schmidt orthogonalization
void CModal::Orth()
{

}

// Integration with the newly G_alpha method
void CG_alpha::G_alpha_Intregration(CLoadCaseData& Load, int i_load)
{
	// Calculate the lumped mass matrix
	Cal_LM();

	// Calculate the lamping matrix C
	int NEQ = M->dim();
	C = new CSkylineMatrix<double>(NEQ);
	C->Generate_C(K, L_M, C_alpha, C_beta);

	// Give the parameter of the NOCH method in Table I of the paper
	double alpha = (2 * rho - 1) / (1 + rho);
	double delta = (3 * rho - 1) / 2.0 / (1 + rho);
	double eta = rho / (1 + rho);
	double epso = eta / (1 + rho);
	double beta = epso / rho;
	double miu = eta;
	double gamma = miu / rho;

	// Give the computing constants
	double m1 = beta * (1.0 - eta) * h * h;
	double m2 = gamma * (1.0 - delta) * h;
	double m3 = 1.0 - alpha;
	double m4 = beta * eta * h * h;
	double m5 = h * (1.0 - alpha);
	double m6 = (epso * (1.0 - alpha) - alpha * beta) * h * h;
	double m7 = (gamma * (1.0 - delta) - beta) * h * h;
	double m8 = (1.0 - delta) * (epso * gamma - beta * miu) * h * h * h;
	double m9 = -beta * eta * h * h;
	double m10 = gamma / beta / h;
	double m11 = (beta - gamma) / beta;
	double m12 = (beta * miu - epso * gamma) * h / beta;
	double m13 = 1.0 / beta / h / h;
	double m14 = -1.0 / beta / h;
	double m15 = -epso / beta;

	// Initialize the dis, vel and acc and them at t - h dis_p, vel_p, acc_p
	double* dis_p = new double[NEQ];
	double* vel_p = new double[NEQ];
	double* acc_p = new double[NEQ];
	double* Force_p = new double[NEQ];
	Recall_message(i_load, dis_p, vel_p, acc_p, Force_p);

	// *** Without C matrix is first tested ***
	double time_load = Load.time_knot;

	// Get the force and acc at time 0	-- Eq(6)
	if (t < eps)
	{
		Force_p = Cur_Force(t, Load);
		for (unsigned int i = 0; i < NEQ; i++){
			acc_p[i] = Force_p[i];
			for (unsigned int j = 0; j < NEQ; j++)
				acc_p[i] -= ((*K)(i + 1, j + 1) * dis[j] - (*C)(i + 1, j + 1) * vel[j]);
			acc_p[i] /= L_M[i];
		}
		int* num_freedom = new int[N_His_Freedom];
		for (unsigned int i = 0; i < N_His_Freedom; i++) {
			num_freedom[i] = NodeList[Freedom_output[2 * i] - 1].bcode[Freedom_output[2 * i + 1] - 1];
		}
		His_Output->OutputHisMessage(t, dis, vel, acc, N_His_Freedom,num_freedom);
	}

	// Genarate the effective stiffness matrix K_e -- Eq(60)
	CSkylineMatrix<double>* K_e = new CSkylineMatrix<double>(NEQ);
	K_e->Generate_Ke(K, C, L_M, m1, m2, m3);

	//// ****for debug****
	//unsigned int* DiagonalAddress = K_e->GetDiagonalAddress();
	//cout << setiosflags(ios::scientific) << setprecision(5);
	//for (int i = 0; i < DiagonalAddress[NEQ] - DiagonalAddress[0]; i++) {
	//	cout << setw(14) << (*K_e)(i);

	//	if ((i + 1) % 6 == 0)
	//	{
	//		cout << endl;
	//	}
	//}
	//cout << endl;
	
	// LDLT the effective stiffness matrix -- Eq(61)
	LDLT(K_e);

	// Count the Tecplot output
	int Tec_Count = 1;

	// Integration
	while (t < time_load)
	{
		t += h;				// t = (i + 1) * h

		double* F = Cur_Force(t, Load);		// F_t, that is, Force
		unsigned int* DiagonalAddress = K->GetDiagonalAddress();

		// Generate the effective force vector -- Eq(62)
		double* F_e = new double[NEQ];
		for (int i = 0; i < NEQ; i++) {
			F_e[i] = m1 * F[i] + m4 * Force_p[i] + L_M[i] * (m3 * dis_p[i] + m5 * vel_p[i] + m6 * acc_p[i]);
			for (int j = 0; j < NEQ; j++) {
				int H;
				if (j >= i)
					H = DiagonalAddress[j + 1] - DiagonalAddress[j];
				else
					H = DiagonalAddress[i + 1] - DiagonalAddress[i];
				if ((j - i - H < 0) && (i - j - H < 0))
				{
					F_e[i] += m9 * (*K)(i + 1, j + 1) * dis_p[j];
					F_e[i] += (*C)(i + 1, j + 1) * (m2 * dis_p[j] + m7 * vel_p[j] + m8 * acc_p[j]);
				}	
			}
		}

		// Back substitution, after calculation F_e will be the dis at current time -- Eq(63)
		BackSubstitution(K_e, F_e);


		// Compute the current vel and acc
		for (unsigned int i = 0; i < NEQ; i++) {
			dis[i] = F_e[i];
			vel[i] = m10 * (dis[i] - dis_p[i]) + m11 * vel_p[i] + m12 * acc_p[i];		//Eq(64)
			acc[i] = m13 * (dis[i] - dis_p[i]) + m14 * vel_p[i] + m15 * acc_p[i];
		}

		// History Output
		int* num_freedom = new int[N_His_Freedom];
		for (unsigned int i = 0; i < N_His_Freedom; i++) {
			num_freedom[i] = NodeList[Freedom_output[2 * i] - 1].bcode[Freedom_output[2 * i + 1] - 1];
		}
		His_Output->OutputHisMessage(t, dis, vel, acc, N_His_Freedom, num_freedom);

		// Tecplot Output
		if (fmod((double)Tec_Count, (double)TecplotOut_Interval) == 0)
			Tecplot_Output->OutputTecplot(t, dis);
		Tec_Count += 1;

		// Store the motion messages of the last step
		Store_message(dis_p, vel_p, acc_p, Force_p);
	}

}

double* CG_alpha::Cur_Force(double time, CLoadCaseData& load)
{
	int NEQ = M->dim();
	clear(Force, NEQ);

	//	Loop over for all concentrated loads in load case LoadCase
	for (unsigned int lnum = 0; lnum < load.nloads; lnum++)
	{
		unsigned int dof = NodeList[load.node[lnum] - 1].bcode[load.dof[lnum] - 1];

		if (dof) // The DOF is activated
			Force[dof - 1] += load.a_[lnum] + load.b_[lnum] * time + load.c_[lnum] * sin(load.d_[lnum] * time + load.e_[lnum]);
	}

	return Force;
}


//! Recall the motion messages of the end of the last load
void CG_alpha::Recall_message(int iload, double* dis_p, double* vel_p, double* acc_p, double* Force_p)
{
	int NEQ = M->dim();
	if (iload == 0)
	{
		dis = new double[NEQ];
		vel = new double[NEQ];
		acc = new double[NEQ];
		Force = new double[NEQ];
		clear(dis, NEQ);
		clear(vel, NEQ);
		clear(acc, NEQ);
		t = 0.0;
	}
	for (unsigned int i = 0; i < NEQ; i++)
	{
		dis_p[i] = dis[i];
		vel_p[i] = vel[i];
		acc_p[i] = acc[i];
		Force_p[i] = Force[i];
	}
}

//! Store the motion messages of the last step
void CG_alpha::Store_message(double* dis_p, double* vel_p, double* acc_p, double* Force_p)
{
	int NEQ = M->dim();
	for (unsigned int i = 0; i < NEQ; i++)
	{
		dis_p[i] = dis[i];
		vel_p[i] = vel[i];
		acc_p[i] = acc[i];
		Force_p[i] = Force[i];
	}
}