/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Bar.h"
#include "Outputter.h"
#include "Clock.h"

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 2) //  Print help message
	{
	    cout << "Usage: stap++ InputFileName\n";
		exit(1);
	}

	string filename(argv[1]);
    size_t found = filename.find_last_of('.');
	string OutFile = filename + ".out";
	CDomain* FEMData = CDomain::Instance();

    // If the input file name is provided with an extension
    if (found != std::string::npos) {
        if (filename.substr(found) == ".dat")
		{
            filename = filename.substr(0, found);
			string InFile = filename + ".dat";

			//  Read data and define the problem domain
			if (!FEMData->ReadData(InFile, OutFile))
			{
				cerr << "*** Error *** Data input failed!" << endl;
				exit(1);
			}
		}
		else if (filename.substr(found) == ".inp")
		{
			filename = filename.substr(0, found);
			string InFile = filename + ".inp";
			OutFile = filename + ".out";
			if (!FEMData->ReadInpData(InFile, OutFile))
			{
				cerr << "*** Error *** Data input failed!" << endl;
				exit(1);
			}
		}
        else {
            // The input file name must has an extension of 'dat'
            cout << "*** Error *** Invalid file extension: "
                 << filename.substr(found+1) << endl;
            exit(1);
        }
    }

    string InFile = filename + ".dat";

	string TecFile = filename + "_tec.dat";

	string vtkFile = filename + ".vtk";

	string HisFile = filename + ".his";

    Clock timer;
    timer.Start();



//  Output the result in tecplot form at the beginning including the head
	COutputter* Tec_Output = COutputter::Tec_Instance(TecFile);
	Tec_Output->OutputTecplot(0);

//  Output head,nodes and elements information to vtkfile
	COutputter* vtk_Output = COutputter::vtk_Instance(vtkFile);
	vtk_Output->OutputVTKHead();
	vtk_Output->OutputVTKNodes();
	vtk_Output->OutputVTKElements();
    
    double time_input = timer.ElapsedTime();

//  Allocate global vectors and matrices, such as the Force, ColumnHeights,
//  DiagonalAddress and StiffnessMatrix, and calculate the column heights
//  and address of diagonal elements
	FEMData->AllocateMatrices();
    
//  Assemble the banded gloabl stiffness matrix
	FEMData->AssembleStiffnessMatrix();

//  Assemble the global mass matrix  1--consistent mass matrix  2--lumped mass matrix
	FEMData->AssembleMassMatrix();
    
    double time_assemble = timer.ElapsedTime();
	COutputter* Output = COutputter::Instance();
	double time_solution = 0.0;

	int Solving_Type = FEMData->GetMODEX();

// ****************************************** //
// ***  For the stastic problem calculate *** //
// ****************************************** //
	if (Solving_Type == 1)
	{ 
	//  Solve the linear equilibrium equations for displacements
		CLDLTSolver* Solver = new CLDLTSolver(FEMData->GetStiffnessMatrix());
    
	//  Perform L*D*L(T) factorization of stiffness matrix
		Solver->LDLT();

	#ifdef _DEBUG_
		Output->PrintStiffnessMatrix();
	#endif
        
	//  Loop over for all load cases
		for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
		{
	//      Assemble righ-hand-side vector (force vector)
			FEMData->AssembleForce(lcase + 1);
            
	//      Reduce right-hand-side force vector and back substitute
			Solver->BackSubstitution(FEMData->GetForce());
            
	#ifdef _DEBUG_
			Output->PrintDisplacement(lcase);
	#endif
            
			Output->OutputNodalDisplacement(lcase);
		}

		time_solution = timer.ElapsedTime();


//  Calculate and output stresses of all elements
	Output->OutputElementStress();

//  Output the result in vtk form
	vtk_Output->OutputVTK();

//  Output the result in tecplot form
	Tec_Output->OutputTecplot(1);

	}


// ******************************** //
// ***  For the modal calculate *** //
// ******************************** //
	else if (Solving_Type == 2)
	{
		//  Modal Analysis
		CModal* Modal_ = new CModal(FEMData->GetStiffnessMatrix());

		//  Lanczos method
		Modal_->Lanczos();

		time_solution = timer.ElapsedTime();

		//  Output

	}


// ******************************************* //
// ***  For the Dynamics problem calculate *** //
// ******************************************* //
	else if (Solving_Type == 3)
	{
		*Output << "S T A R T: D Y N A M I C S  A N A L Y S I S" << endl << endl;
		// Dynamics analysis
		CG_alpha* G_alpha_ = new CG_alpha(FEMData->GetStiffnessMatrix(), FEMData->GetMassMatrix());

		CLoadCaseData* Loads = FEMData->GetLoadCases();

		G_alpha_->Obtain_NodeList(FEMData->GetNodeList());

		G_alpha_->Obtain_Dyn_Para(FEMData->GetDynPara());

		//  Output the history result of certain freedom
		COutputter* His_Output = COutputter::His_Instance(HisFile);
		G_alpha_->Obtain_HisOutput(His_Output, FEMData->GetNumHisFreedom(), FEMData->GetMessHisFreedom());

		// Tecplot Output
		G_alpha_->Obtain_TecOutput(Tec_Output);

		// Paraview Output
		G_alpha_->Obtain_VTKOutput(vtk_Output);

		for (unsigned int i = 0; i < FEMData->GetNLCASE(); i++)
		{
			*Output << "	Begin the Load case		" << i + 1 << endl;
			// Integration with the newly G_alpha method
			G_alpha_->G_alpha_Intregration(Loads[i], i, filename);
			*Output << "	Finish the Load case		" << i + 1 << endl;
		}



		time_solution = timer.ElapsedTime();

	}
// ******************************************* //

	else
	{
		*Output << "E R R: No such solving type!";
	}

    double time_stress = timer.ElapsedTime();
    
    timer.Stop();
    
    *Output << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
            << "     TIME FOR INPUT PHASE = " << time_input << endl
            << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl
            << "     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS = " << time_solution - time_assemble << endl << endl
            << "     T O T A L   S O L U T I O N   T I M E = " << time_stress << endl;

	return 0;
}
