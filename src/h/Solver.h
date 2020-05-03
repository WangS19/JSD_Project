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

#include "SkylineMatrix.h"
#include "LoadCaseData.h"
#include "Node.h"
#include "Outputter.h"
#include <vector>
#include <algorithm>

//!	Clear an array
template <class type> void clear(type* a, unsigned int N);

//!	Base class for a solver
/*	New solver should be derived from this base class, and match the storage scheme
	of the global stiffness matrix employed in Domain class. */
class CSolver
{
protected:

	CSkylineMatrix<double>* K;

	CSkylineMatrix<double>* M;

	CSkylineMatrix<double>* C;

	// The lumped mass matrix
	double* L_M;

public:

//! Constructor
	CSolver(CSkylineMatrix<double>* K);

//! Overloading the constructing function
	CSolver(CSkylineMatrix<double>* K, CSkylineMatrix<double>* M);		// Overloading the constructing function

//! Calculate the lumped mass matrix
	void Cal_LM();

//! Perform L*D*L(T) factorization of any of the matrix with skyline form
	void LDLT(CSkylineMatrix<double>* KK);

//! Back substitute
	void BackSubstitution(CSkylineMatrix<double>* KK, double* Force);
    
};

//!	LDLT solver: A in core solver using skyline storage  and column reduction scheme
class CLDLTSolver : public CSolver
{
public:

//!	Constructor
	CLDLTSolver(CSkylineMatrix<double>* K) : CSolver(K) {};

//!	Perform L*D*L(T) factorization of the stiffness matrix
	void LDLT();

//!	Reduce right-hand-side load vector and back substitute
	void BackSubstitution(double* Force); 
};

//! Modal solver: modal solver with Lanczos method
class CModal : public CSolver
{
private:

	//! The tolerance for Jacobi
	double Tol_J; //= 1.0e-9;

	//! The tolerance for Lanczos
	double Tol_L; //= 1.0e-9;

	//! Number of required eigenvalues *(Input)
	int Nroot; //= 2;

	//! Maximum number of the restart, set 5 usually
	int N_ite_max; //= 5;

	//! Number of itegration vectors used
	//! Usually set to be min(2*Nroot, Nroot+8), but less than the freedom of system
	int N_iv ;

	//! The eigenvalue vector
	double* Eig_value;

	//! The eigenvector matrix
	std::vector<double>* Eig_vector;



public:

	//! Constructor
	CModal(CSkylineMatrix<double>* K) : CSolver(K),Tol_J(1.0e-9),Tol_L(1.0e-9),Nroot(2),N_ite_max(5) {};

	//! To solve the generalized eigenproblem with the Lanczos method
	void Lanczos();

	//! To solve the standard eigenproblem with the Jacobi method
	void Jacobi();

	//! Use the Gram-Schmidt orthogonalization
	void Orth();

};

//! Modal solver: modal solver with Lanczos method
class CG_alpha : public CSolver
{
private:
	//! spectral radius

	double rho = 0.5;

	//! Damping coefficient *******
	double C_alpha;
	double C_beta;


	//! The current motion message
	double* dis;
	double* vel;
	double* acc;

	//! The current time
	double t;

	//! The current force
	double* Force;

	//! The time step	*******

	double h;

	//! The minimal quantity
	const double eps = 1e-16;

	//!	List of all nodes in the domain
	CNode* NodeList;

	//! History Output
	COutputter* His_Output;

	//! Tecplot Output
	COutputter* Tecplot_Output;

	int TecplotOut_Interval = 20;


	//! History output which freedom
	int N_His_Freedom;
	int* Freedom_output;

public:

	//! Constructor

	CG_alpha(CSkylineMatrix<double>* K, CSkylineMatrix<double>* M) : CSolver(K,M) {};


	//! Integration
	void G_alpha_Intregration(CLoadCaseData& Load, int i_load);

	//! Recall the motion messages of the end of the last load
	void Recall_message(int iload, double* dis, double* vel, double* acc, double* Force);

	//! Store the motion messages of the last step
	void Store_message(double* dis, double* vel, double* acc, double* Force);

	//! Calculate the external force at the current time
	double* Cur_Force(double t, CLoadCaseData& load);

	//! Obtain the node list message
	void Obtain_NodeList(CNode* N_list) { NodeList = N_list; }

	//! Obtain the histoty output file
	void Obtain_HisOutput(COutputter* H_Output, int N_History_Freedom, int* Freedom_message_output) {
		His_Output = H_Output; 
		N_His_Freedom = N_History_Freedom;
		Freedom_output = Freedom_message_output;
	}

	//! Obtain the Tecplot output file
	void Obtain_TecOutput(COutputter* Tec_Output) { Tecplot_Output = Tec_Output; }

	//! Obtain the dynamics parameters
	void Obtain_Dyn_Para(double* Dyn_Para) {
		h = Dyn_Para[0];
		C_alpha = Dyn_Para[1];
		C_beta = Dyn_Para[2];
	}

};