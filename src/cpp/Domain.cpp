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
#include "Material.h"

using namespace std;

//	Clear an array
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;
	
	NUMEG = 0;
	EleGrpList = nullptr;
	
	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;
	
	NEQ = 0;

	Force = nullptr;
	StiffnessMatrix = nullptr;

}

//	Desconstructor
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::Instance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

// Read domain data from the inp file
bool CDomain::ReadInpData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::Instance(OutFile);

	//Read the headline

	Input.seekg(0);
	Input.getline(Title,256);
	string line1;
	string line2;
	getline(Input,line1);
	getline(Input,line2);
	strcpy(Title,(line1+"\n"+line2).c_str() );
	Output->OutputHeading();


	// Read NUMNP,NUMEG,NLCASE,MODEX
	int flag=0;
	int *pflag= & flag;
	
	streampos pnode = seek("*Node",0,pflag);
	if(flag)
		cerr << "*** ERROR *** No Nodes Information!" << endl;
	streampos pelement = seek("*Element",pnode,pflag);
	if(flag)
		cerr << "*** ERROR *** No Element Information!" << endl;
	NUMNP = countcode( " ", pnode, pelement);

	streampos psolidtitle = seek( "** Section:", 0, pflag);
	streampos pendpart = seek( "*End Part",psolidtitle, pflag);
	NUMEG = countcode( "*Solid Section", psolidtitle , pendpart);
	streampos psteptitle = seek("** STEP",0,pflag);
	if(flag)
		cerr << "*** ERROR *** No Load Step Information!" << endl;
	Input.seekg(-1, ios :: end);
	streampos pend = Input.tellg();
	NLCASE = countcode("*Step", psteptitle, pend);
	
	streampos pstatic = seek("*Static",0,pflag);
	if(!flag)
		MODEX = 1;
	else
	{
		streampos pdynamic = seek("*Dynamic",0 , pflag);
		if(!flag)
			MODEX = 3;
		else
			MODEX = 2;
//			cerr<< "*** ERROR *** No Loads Information!" << endl;
	}

	//Read Nodal Point Data
	if (ReadInpNodalPoints())
        Output->OutputNodeInfo();
    else
        return false;

	//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

	//  Read the dynamics parameters

	if (MODEX == 3)
	{
		Dyna_para[1] = 0.01;
		Dyna_para[2] = 0.01;
		streampos pdynamic = seek( "*Dynamic",0 , pflag);
		Input.seekg( pdynamic, ios_base :: beg);
		size_t comma = 0;
		string line;
		getline(Input,line);
		getline(Input,line);
		comma = line.find(',',0);
		Dyna_para[0] = atof( line.substr(0, comma).c_str());
	}

	//	Read load data
	if (ReadInpLoadCases())
        Output->OutputLoadInfo();
    else
        return false;

//	Read element data
	if (ReadInpElements())
        Output->OutputElementInfo();
    else
        return false;

//	Read history output message
	if (MODEX == 3)
		if (!ReadInpHisMessage())
			return false;

//  Read Animation output message

	return true;
}
streampos CDomain :: seekset(string target, string set)
{
	//seek set-?
	int flag = 0;
	int* pflag = &flag;
	streampos passemblytitle = seek( "** ASSEMBLY",0,pflag);
	streampos pmaterialtitle = seek( "** MATERIALS",0,pflag);
	Input.seekg(passemblytitle, ios_base :: beg);

	string line;
	int setlength=set.length();
	char *buff = new char[setlength];
	memset(buff,0,setlength+1);
	Input.read(buff,setlength);
	char *content = new char [setlength];
	strcpy(content,set.c_str());
	while(Input.read(buff,setlength),Input.tellg() < pmaterialtitle)
	{
		Input.seekg( -setlength , ios :: cur );
		getline(Input,line);
		if(!strcmp(content,buff))
		{
			size_t commaa = 0;
			size_t commab = 0;
			commaa=line.find(',',0);
			commab=line.find(',',commaa+1);
			if (commab > line.length() )
				commab =line.length();
			string tmp = line.substr(commab-target.length(),target.length()).c_str();
			if (!(stricmp(tmp.c_str(), target.c_str())))                        // find "Set-?"
				break;
		}
	}
	return Input.tellg();
}

streampos CDomain::seek(string code, streampos start , int *pflag)
{
	Input.seekg(start, ios_base:: beg);
	string line;
	*pflag = 0;
	int length = code.length();
	char *content = new char [length];
	strcpy(content,code.c_str());
	char *tmp = new char [length];
	memset(tmp,0,length+1);

	while(Input.read(tmp,length),!Input.eof())
	{
		Input.seekg(-length,ios :: cur);
		if(!strcmp(tmp,content))
			break;
		getline(Input,line);
	}
	if( Input.eof())
	{
		Input.clear();  //clear eof!!!
		Input.seekg(0, ios_base::beg);
		*pflag = 1;
		return 0;
	}
	streampos p=Input.tellg();
	return p;
}

int CDomain::countcode(string code, streampos start, streampos end)
{
	Input.seekg(start, ios :: beg);
	streampos p=Input.tellg();
	string line;
	int length = code.length();
	char *content = new char [length];
	strcpy(content,code.c_str());
	char *tmp = new char [length];
	memset(tmp,0,length+1);
	int count = 0;
	
	while(p < end)
	{
		Input.read(tmp,length);
		Input.seekg(-length,ios :: cur);
		if(!strcmp(tmp,content))
			count ++;
		getline(Input,line);
		p=Input.tellg();
	}
	return count;
}

//not finished,read BC
bool CDomain::ReadInpNodalPoints()
{
	//	Read nodal point data lines
	NodeList = new CNode[NUMNP];
	string line;
	int flag = 0;
	int *pflag = & flag;

	streampos pnode = seek("*Node",0,pflag);
	if(flag)
		cerr << "*** ERROR *** No Nodes Information!" << endl;
	Input.seekg(pnode);
	getline(Input,line);

	char *c= new char[1];
	memset(c,0,2);
	int np = 0;

	// Read nodal coordinate information
	Input.read(c,1);
	while(c[0]!='*')//if c[0]=='*',Nodes info has been read.
	{
		Input.seekg(-1, ios :: cur);
		getline(Input,line);
		size_t comma=0;
		size_t comma2=0;
		int count=0;
		double data[3]={0};//innitial coordinate to 0,0,0;


		comma=line.find(',',0);//find the place of first comma
		NodeList[np].NodeNumber = atoi(line.substr(0,comma).c_str());//get string from 0 to comma. change format into int.
		if(NodeList[np].NodeNumber == 0)
		{
			cerr << "*** ERROR *** Nodes Information ERROR!" << endl;
			return false;
		}
		while( comma < line.size())//until the end of the line
		{
			comma2 = line.find(',',comma+1);//find next comma
			count++;// count to judge dimension, 1 for 1D, 2 for 2D, 3 for 3D
			if (comma2 >= line.size())//no more commas
				comma2 = line.size();//set comma2 to the end of line
			data[count-1] = atof(line.substr(comma+1,comma2-comma-1).c_str());
			comma= comma2;
		}
		Input.read(c,1);
		switch(count)
		{
			case 1:
				NodeList[np].XYZ[0] = data[0];
				NodeList[np].XYZ[1] = 0;
				NodeList[np].XYZ[2] = 0;
				break;
			case 2:
				NodeList[np].XYZ[0] = data[0];
				NodeList[np].XYZ[1] = data[1];
				NodeList[np].XYZ[2] = 0;
				break;
			case 3:
				NodeList[np].XYZ[0] = data[0];
				NodeList[np].XYZ[1] = data[1];
				NodeList[np].XYZ[2] = data[2];
				break;
		}
		np++;
	}

	//Initial BC
	for (int i = 0 ; i < NUMNP ; i ++)
	{
		for ( int j=0 ; j < 2 ; j ++)
			NodeList[i].bcode[j] = 0;
		NodeList[i].bcode[2] = 1;
	}

	//Read Boundary Condition
	streampos pboundary = seek("*Boundary",0,pflag);
	streampos ploadtitle = seek("** LOADS",0,pflag);
	Input.seekg( pboundary, ios_base::beg);
	getline(Input,line);
	memset(c,0,2);
	Input.read(c,1);
	while(Input.tellg() < ploadtitle)
	{
		Input.seekg(-1,ios :: cur);
		getline(Input,line);
		if(c[0]!='*')
		{
			size_t comma=0;
			size_t comma2=0;
			int count=0;
			int data[3]={0};//innitial coordinate to 0,0,0;
			string set;
			

			comma=line.find(',',0);//find the place of first comma
			set = line.substr(0,comma).c_str();//get string from 0 to comma. set="Set-?"
			
			while( comma < line.size())//until the end of the line
			{
				comma2 = line.find(',',comma+1);//find next comma
				count++;// count to judge dimension, 1 for 1BC, 2 for no1-no2, 3 for no1-no2 per no3
				if (comma2 >= line.size())//no more commas
					comma2 = line.size();//set comma2 to the end of line
				data[count-1] = atoi(line.substr(comma+1,comma2-comma-1).c_str());
				comma= comma2;
			}
			
			//seek "Set-?" information between "*End Instance" & "*End Assembly"

			streampos pnow = Input.tellg();
			streampos pend_instance = seek( "*End Instance",0,pflag);
			streampos pend_assembly = seek( "*End Assembly",0,pflag);
			Input.seekg(pend_instance, ios_base :: beg);
			//Search "*Nset"
			string code = "*Nset";
			int codelength=code.length();
			char *buff = new char[codelength];
			memset(buff,0,codelength+1);
			Input.read(buff,codelength);
			char *content = new char [codelength];
			strcpy(content,code.c_str());
			while(Input.read(buff,codelength),Input.tellg() < pend_assembly)
			{
				Input.seekg( -codelength , ios :: cur );
				getline(Input,line);
				if(!strcmp(content,buff))
				{
					size_t commaa = 0;
					size_t commab = 0;
					commaa=line.find(',',0);
					commab=line.find(',',commaa+1);
					string tmp = line.substr(commaa+7,5).c_str();
					if (tmp == set)                        // find "Set-?"
					{
						while(Input.read(buff,codelength),Input.tellg() < pend_assembly)  //Judging read info or not
						{
							Input.seekg( -codelength , ios :: cur );
							getline(Input,line);
							if(!strcmp(content,buff))  // searching reach next "*Nset",break
								break;

							// read info, *.inp allows not more than 16 data in one line.

							int NBC[16] = {0};
							commaa=line.find(',',0);
							NBC[0] = atoi(line.substr(0,commaa).c_str());
							int countnp = 1;
							while( commaa < line.size()-1)
							{
								countnp++;
								commab=line.find(',',commaa+1);
								if (commab >= line.size())//no more commas
									commab = line.size();//set commab to the end of line
								NBC[countnp-1] = atoi(line.substr(commaa+1,commab-commaa-1).c_str());
								commaa = commab;
							}

							switch(count)
							{
							case 1:    // 1 node BC
								if(data[0] <= 3)
								{
									int count1 = 0;
									while( count1 < countnp)
									{
										NodeList[NBC[count1]].bcode[data[0]] = 1;
										count1 ++;
									}
								}
								break;
							case 2:    // from No1 to No2 BC
								if(data[0] <= 3)
								{
									int count1 = 0;
									while( count1 < countnp)
									{
										for( int k = data[0] ; k <= data[1]; k++)
											NodeList[NBC[count1]-1].bcode[k-1] = 1;
										count1 ++;										
									}
								}
								break;
							case 3:    // from No1 to No2 per No3 BC
								if(data[0] <= 3)
								{
									int count1 = 0;
									while( count1 < countnp)
									{
										for( int k = data[0] ; k <= data[1]; k += data[2])
											NodeList[NBC[count1]].bcode[k] = 1;
										count1 ++;										
									}
								}
								break;
							}
						}
					}
				}
			}
			
			Input.seekg(pnow, ios_base :: beg);				
		}
		Input.read(c,1);
	}	
	return true;
}

bool CDomain::ReadInpLoadCases()
{
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases
	
	string line;
	int flag = 0;
	int *pflag = & flag;
	streampos pnow = 0;
	streampos pstep = seek("*Step",pnow,pflag);
	streampos pendstep = seek("*End Step",pstep,pflag);
	unsigned int NCL = countcode("*Cload",pstep,pendstep);
	unsigned int NL = 0;
	streampos pend_assembly = seek( "*End Assembly",0,pflag);

	switch(MODEX)
	{
	case 1://Static

		for (unsigned lcase = 0; lcase < NLCASE ; lcase ++)
		{

			NCL = countcode("*Cload",pstep,pendstep);
			NL = 0;

			for (unsigned int lnum = 0; lnum < NCL ; lnum++)
			{
				pnow = seek("*Step", 0, pflag);
				for (unsigned int nstep = 0; nstep < lcase; nstep ++)
				{
					getline(Input,line);
					pnow = Input.tellg();
					pnow = seek("*Step", pnow, pflag);
				}
				Input.seekg(pnow, ios :: beg);
				for (unsigned int cload = 0; cload < lnum + 1; cload ++)
				{
					getline(Input,line);
					pnow = Input.tellg();
					pnow = seek("*Cload",pnow ,pflag);
				}
				Input.seekg(pnow, ios :: beg);
				getline(Input,line);
				size_t comma = 0;
				size_t comma1 = 0;
				getline(Input,line);
				comma = line.find(',',0);
				comma1 = line.find(',',comma+1);
				string set = line.substr(0,comma).c_str();    // read set-?
				pnow = Input.tellg();

				streampos pthisset = seekset( set, "*Nset");
				Input.seekg(pthisset, ios_base :: beg);
				char *buff = new char [5];
				memset(buff,0,6);
				char *content = new char [5];
				strcpy(content,"*Nset");
				while(Input.read(buff,5))
				{
					Input.seekg(-5, ios :: cur);

					if(!strcmp(content,buff))
						break;
					getline(Input,line);
					if(Input.tellg() > pend_assembly)
						break;
					size_t commaa = 0;
					while(commaa < line.length()-1)
					{
						commaa = line.find(',',commaa+1);
						NL++;
					}
				}

			}
			Input.seekg(pnow, ios_base :: beg);
			LoadCases[lcase].Allocate(NL);
			unsigned fnum = 0;

			pnow = seek("*Step", 0, pflag);
				for (unsigned int nstep = 0; nstep < lcase; nstep ++)
				{
					getline(Input,line);
					pnow = Input.tellg();
					pnow = seek("*Step", pnow, pflag);
				}
			pnow = Input.tellg();
			for (unsigned int lnum = 0; lnum < NCL ; lnum++)
			{
				Input.seekg(pstep, ios_base :: beg );
				pnow = seek("*Cload",pnow ,pflag);
				getline(Input,line);
				size_t comma = 0;
				size_t comma1 = 0;
				getline(Input,line);
				comma = line.find(',',0);
				comma1 = line.find(',',comma+1);
				string set = line.substr(0,comma).c_str();    // read set-?
				pnow = Input.tellg();

				int dof = atoi(line.substr(comma+1,comma1-comma-1).c_str());
				double load = atof(line.substr(comma1+1,line.size()-comma1-2).c_str());

				streampos pthisset = seekset( set, "*Nset");
				Input.seekg( pthisset, ios_base :: beg);
				char *c = new char[1];
				memset(c,0,2);
				Input.read(c,1);
				while (strcmp(c,"*"))
				{
					Input.seekg(-1, ios :: cur);
					getline(Input,line);
					size_t commaa = 0;
					size_t commab = 0;
					commaa = line.find(',',0);
					LoadCases[lcase].node[fnum] = atoi(line.substr(0,comma).c_str());
					LoadCases[lcase].dof[fnum] = dof;
					LoadCases[lcase].load[fnum] = load;
					while ( commaa < line.length())
					{
						fnum++;
						commab = line.find(',',commaa+1);
						if ( !(commaa == line.length() - 1) )
						{
							commab = line.find(',',commaa+1);
							if (commab > line.length())
								commab = line.length();
							LoadCases[lcase].node[fnum]  = atoi(line.substr(commaa+1,commab-commaa-1).c_str());
							LoadCases[lcase].dof[fnum] = dof;
							LoadCases[lcase].load[fnum] = load;
							commaa = commab;
						}
						else
							break;
					}
					Input.read(c,1);
				}
			}

		}
		break;
	case 2:
		break;
	case 3: //dynamic
		for (unsigned lcase = 0; lcase < NLCASE ; lcase ++)
		{

			NCL = countcode("*Cload",pstep,pendstep);
			NL = 0;

			for (unsigned int lnum = 0; lnum < NCL ; lnum++)
			{
				pnow = seek("*Step", 0, pflag);
				for (unsigned int nstep = 0; nstep < lcase; nstep ++)
				{
					getline(Input,line);
					pnow = Input.tellg();
					pnow = seek("*Step", pnow, pflag);
				}
				for (unsigned int cload = 0; cload < lnum + 1; cload ++)
				{
					getline(Input,line);
					pnow = Input.tellg();
					pnow = seek("*Cload",pnow ,pflag);
				}
				getline(Input,line);
				size_t comma = 0;
				size_t comma1 = 0;				
				getline(Input,line);
				comma = line.find(',',0);
				comma1 = line.find(',',comma+1);
				string set = line.substr(0,comma).c_str();    // read set-?
				pnow = Input.tellg();

				streampos pthisset = seekset( set, "*Nset");
				Input.seekg(pthisset, ios_base :: beg);
				char *buff = new char [1];
				memset(buff,0,2);
//				char *content = new char [5];
//				strcpy(content,"*Nset");
				while(Input.read(buff,1))
				{
					Input.seekg(-1, ios :: cur);

					if(!strcmp("*",buff))
						break;
					getline(Input,line);
					if(Input.tellg() > pend_assembly)
						break;
					size_t commaa = 0;
					while(commaa < line.length()-1)
					{
						commaa = line.find(',',commaa+1);
						NL++;
					}
				}

			}
			Input.seekg(pnow, ios_base :: beg);
			LoadCases[lcase].Allocate(NL);
			unsigned int fnum = 0;

			pnow = seek("*Step", 0, pflag);
				for (unsigned int nstep = 0; nstep < lcase; nstep ++)
				{
					getline(Input,line);
					pnow = Input.tellg();
					pnow = seek("*Step", pnow, pflag);
				}
			pnow = Input.tellg();
			for (unsigned int lnum = 0; lnum < NCL ; lnum++)
			{
				Input.seekg(pstep, ios_base :: beg );
				pnow = seek("*Cload",pnow ,pflag);
				getline(Input,line);
				size_t comma = 0;
				size_t comma1 = 0;
				size_t pequal = 0;
				pequal = line.find('=',0);
				string amplitude = line.substr(pequal+1,line.length());
				getline(Input,line);
				comma = line.find(',',0);
				comma1 = line.find(',',comma+1);
				string set = line.substr(0,comma).c_str();    // read set-?
				pnow = Input.tellg();

				int dof = atoi(line.substr(comma+1,comma1-comma-1).c_str());
				double load = atof(line.substr(comma1+1,line.size()-comma1-2).c_str());
				double a_ = 0;
				double b_ = 0;
				double c_ = 0;
				double d_ = 0;
				double e_ = 0;

				Dyna_para[1] = 0.01;
				Dyna_para[2] = 0.01;
				streampos pdynamic = seek( "*Dynamic",0 , pflag);
				Input.seekg( pdynamic, ios_base :: beg);
				size_t comma_1 = 0;
				size_t comma_2 = 0;
				getline(Input,line);
				getline(Input,line);
				comma_1 = line.find(',',0);
				comma_2 = line.find(',',comma_1+1);
				LoadCases[lcase].time_knot = atof(line.substr(comma_1+1,comma_2-comma_1-1).c_str());

				streampos pamp = seekset( amplitude,"*Amplitude");
				Input.seekg(pamp, ios_base :: beg);
				char *cbuff = new char [1];
				memset(cbuff, 0 , 2);
				while(Input.read(cbuff,1))
				{
					Input.seekg(-1, ios :: cur);
					getline(Input,line);
					if(!strcmp(cbuff,"*"))
						break;
					size_t commaa = 0;
					size_t commab = 0;
					commaa = line.find(',',0);
					commab = line.find(',',commaa+1);
					double t1 = atof(line.substr(0,commaa).c_str());
					double f1 = atof(line.substr(commaa+1,commab-commaa-1).c_str());
					commaa = commab;
					commab = line.find(',',commaa+1);
					double t2 = atof(line.substr(commaa+1,commab-commaa-1).c_str());
					commaa = commab;
					commab = line.find(',',commaa+1);
					if(commab > line.length())
						commab = line.length();
					double f2 = atof(line.substr(commaa+1,commab-commaa-1).c_str());
					b_ = (f1 - f2)/ (t1 - t2);
					a_ = f1 - b_ * t1;
				}

				streampos pthisset = seekset( set, "*Nset");
				Input.seekg( pthisset, ios_base :: beg);
				char *cbuff1 = new char[1];
				memset(cbuff1,0,2);
				Input.read(cbuff1,1);
				while (strcmp(cbuff1,"*"))
				{
					Input.seekg(-1, ios :: cur);
					size_t commaa = 0;
					streampos ptemp = Input.tellg();
					getline(Input,line);
					unsigned int length = line.length();
					while(commaa < length)
					{
						commaa = line.find(',',0);
						Input.seekg(ptemp, ios_base::beg);
						getline(Input,line,',');
						ptemp = Input.tellg();
						LoadCases[lcase].node[fnum] = atoi(line.c_str());
						LoadCases[lcase].dof[fnum] = dof;
						LoadCases[lcase].load[fnum] = load;
						LoadCases[lcase].a_[fnum] = a_;
						LoadCases[lcase].b_[fnum] = b_;
						LoadCases[lcase].c_[fnum] = c_;
						LoadCases[lcase].d_[fnum] = d_;
						LoadCases[lcase].e_[fnum] = e_;
						getline(Input,line);
						length = line.length();
						fnum ++;
					}
					LoadCases[lcase].node[fnum] = atoi(line.c_str());
					LoadCases[lcase].dof[fnum] = dof;
					LoadCases[lcase].load[fnum] = load;
					LoadCases[lcase].a_[fnum] = a_;
					LoadCases[lcase].b_[fnum] = b_;
					LoadCases[lcase].c_[fnum] = c_;
					LoadCases[lcase].d_[fnum] = d_;
					LoadCases[lcase].e_[fnum] = e_;
					fnum ++;
					Input.read(cbuff1,1);
				}
			}

		}
		break;
		
	}

	return true;
}

bool CDomain::ReadInpElements()
{
	EleGrpList = new CElementGroup[NUMEG];

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		streampos pnow = 0;
		int flag = 0;
		int * pflag = &flag;
		for (unsigned int i = 0; i < EleGrp+1 ; i++)
			pnow = seek("*Solid Section", pnow , pflag);
		string line;
		getline(Input,line);
			size_t comma = 0;
		size_t pequal = 0;
		comma = line.find(',', 0);
		comma = line.find(',', comma + 1);
		pequal = line.find('=', 0);
		string eset = line.substr(pequal + 1, comma - pequal - 1);
		pequal = line.find('=', pequal + 1);
		string mset = line.substr(pequal + 1, line.length() - pequal - 1);
		string etype = " ";
		Input.seekg( seek("*Element",0,pflag), ios_base :: beg);
		size_t posequal = 0;
		getline(Input,line);
		posequal = line.find('=',0);
		etype = line.substr(posequal + 1, line.length() - posequal -1);
		unsigned int NUME = 0;

	//seek position of Material
		
		streampos pMaterialtitle = seek ("** MATERIALS", 0 , pflag);
		streampos pMaterial = seek ("*Material", pMaterialtitle , pflag);
		streampos psteptitle = seek ("** STEP", 0 , pflag);
		Input.seekg( pMaterial , ios_base :: beg);
		*pflag = 0;

		getline(Input,line);
		pequal = line.find('=' , 0);
		unsigned int m_set = 1;
		string tmp = line.substr( pequal+1, line.length() - pequal - 1);
		while (!(tmp == mset))
		{
			m_set ++;
			streampos pMaterial = 0;
			pMaterial = seek ("*Material", pMaterial , pflag);
			Input.seekg( pMaterial , ios_base :: beg);
			getline(Input,line);
			pequal = line.find('=' , 0);
			string tmp = line.substr( pequal+1, line.length() - pequal - 1);
		}

		unsigned int enodes = 0;

		if ( etype == "T2D2")
			enodes = 2;
		else if (etype == "CPE4")
			enodes = 4;
		else if (etype == "CAX8R")
			enodes = 8;
		

		

	//seek ElementType_
		if ( EleGrp == 0)
		{
			streampos pelement = seek ("*Element" , 0 ,pflag);
			Input.seekg( pelement, ios_base :: beg);
			getline(Input, line);
			pequal = line.find('=', 0);
			etype = line.substr(pequal + 1, line.length() - pequal - 1);

			// count NUME
			Input.seekg( pelement, ios_base :: beg);
			getline(Input,line);
			char *c = new char [1];
			memset(c,0,2);
			Input.read(c,1);
			while( stricmp(c,"*"))
			{
				Input.seekg( -1 , ios :: cur);
				NUME ++;
				getline(Input,line);
				Input.read(c,1);
			}
		}
		else
		{
		}
		unsigned int ** Elements = new unsigned int *[NUME];
		for (int i = 0; i < NUME; i++)
			Elements[i] = new unsigned [enodes];

		//Find elset 

		streampos pelement = seek("*Elset",0,pflag); // find elset position
		streampos pele = seek( "*Element",0 , pflag);
		streampos pend = seek( "** ASSEMBLY",0,pflag);
		Input.seekg(pelement, ios_base :: beg);
		size_t pequaltmp = 0;
		while( Input.tellg() < pend)
		{
			getline(Input,line);
			pequaltmp = line.find('=',0);
			size_t comma3 = 0;
			comma3 = line.find(',',comma3+1);
			comma3 = line.find(',',comma3+1);
			if (comma3 > line.length())
				comma3 = line.length();
			if(!(stricmp(eset.c_str(),line.substr(pequaltmp+1,comma3-pequaltmp-1).c_str())))
			{
				pelement = Input.tellg();
				break;
			}
		}

		Input.seekg(pelement, ios_base :: beg);

		//Judge data format 1,16,1/1,2,3,4
		getline(Input,line);
		int length = line.length();
		Input.seekg( - length , ios :: cur );
		size_t comma_tmp = 0;
		comma_tmp = line.find(',',0);
		unsigned int countcomma = 0;
		unsigned int elmode = 0;
		
		while (comma_tmp < line.length())
		{
			countcomma++;
			comma_tmp = line.find(',',comma_tmp + 1);
		}// count comma
	
		if (countcomma == 2 | countcomma == 3)
		{
			getline(Input,line,',');
			getline(Input,line,',');
			unsigned int ntmp1 = 0;
			ntmp1 = atoi(line.c_str());
			if( countcomma == 2)
				getline(Input,line);
			if( countcomma == 3)
				getline(Input,line,',');
			unsigned int ntmp2 = 0;
			ntmp2 = atoi(line.c_str());
			if (ntmp2 < ntmp1)
				elmode = 1;//type 1,16,1 
		}
		
		//Read element num
		unsigned int * elnumset = new unsigned int [ NUME ];
		for (unsigned int i = 0; i < NUME ; i ++)
			elnumset[i] = 0;
		unsigned int elnum = 0;
		char *csign = new char [1];
		memset(csign, 0, 2);
		Input.seekg(pelement, ios_base :: beg);
		Input.read(csign,1);
		
		switch ( elmode )
		{
		case 0:
			//Judge whether data or not
			
			
			while((stricmp(csign,"*")))
			{
				Input.seekg(-1, ios :: cur);
				streampos ptmp = Input.tellg();
				// count comma
				countcomma = 0;
				comma_tmp = line.find(',',0);
				while (comma_tmp < line.length())
				{
					countcomma++;
					comma_tmp = line.find(',',comma_tmp + 1);
				}// count comma finished
				//back to read data
				Input.seekg(ptmp, ios_base :: beg);
				unsigned int elnumtmp = 0;
				for(unsigned int i = 0; i < countcomma; i++)
				{
					getline(Input,line,',');
					elnumtmp = atoi(line.c_str());
					elnumset[ elnum ] = elnumtmp;
					elnum ++;
				}
				getline (Input,line);
				elnumtmp = atoi(line.c_str());
				if (!elnumtmp)// data line is not ended by ','
				{
					elnumset[ elnum ] = elnumtmp;
					elnum ++;
				}

				Input.read(csign,1);
			}
			break;//type 1,2,3,4
		case 1:
			unsigned int tmp1,tmp2,tmp3;
			getline(Input,line,',');
			tmp1 = atoi(line.c_str());
			getline(Input,line,',');
			tmp2 = atoi(line.c_str());
			getline(Input,line);
			tmp3 = atoi(line.c_str());
			for ( unsigned int i = tmp1; i <= tmp2; i += tmp3)
			{
				elnumset[ elnum ] = i;
				elnum ++;
			}
			break;//type 1,16,1
		}
		
		//seek node num of each element, write into Elements[][]
		streampos pseekelement = seek ("*Element",0,pflag);
		
		for (unsigned int i = 0; i < NUME; i++)
		{
			Input.seekg( pseekelement , ios_base :: beg);
			getline(Input,line);
			while(1)
			{
				streampos ptmp = Input.tellg();
				getline(Input,line,',');
				unsigned int tmp1 = atoi(line.c_str());
				if (tmp1 == elnumset [ i ])
				{
					for (unsigned int j = 0; j < enodes - 1; j ++)
					{
						getline(Input,line,',');
						tmp1 = atoi(line.c_str());
						Elements[ elnumset [ i ] - 1 ] [ j ] = tmp1;
					}
					getline(Input,line);
					tmp1 = atoi(line.c_str());
					Elements[ elnumset [ i ] - 1 ] [ enodes - 1 ] = tmp1;
					break;
				}
				getline(Input,line);
			}
		}

        if (!EleGrpList[EleGrp].ReadInp(Input,etype,NUME,pMaterial,m_set,Elements))
            return false;
		for (int i = 0; i < NUME; i++)
			delete[]Elements[i];
		delete[]Elements;
	}
    return true;
}

bool CDomain::ReadInpHisMessage()
{
	Num_His_Output = 1;
	His_freedom = new int[Num_His_Output * 2];
	for (int i = 0; i < Num_His_Output * 2; i++) {
		His_freedom[i] = 1;
	}
	return true;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::Instance(OutFile);

//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();

//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data
	if (ReadNodalPoints())
        Output->OutputNodeInfo();
    else
        return false;

//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

//  Read the dynamics parameters
	if (MODEX == 3)
		Input >> Dyna_para[0] >> Dyna_para[1] >> Dyna_para[2];

//	Read load data
	if (ReadLoadCases())
        Output->OutputLoadInfo();
    else
        return false;

//	Read element data
	if (ReadElements())
        Output->OutputElementInfo();
    else
        return false;

//	Read history output message
	if (MODEX == 3)
		if (!ReadHisMessage())
			return false;

//  Read Animation output message

	return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{

//	Read nodal point data lines
	NodeList = new CNode[NUMNP];

//	Loop over for all nodal points
	for (unsigned int np = 0; np < NUMNP; np++)
		if (!NodeList[np].Read(Input, np))
			return false;

	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof]) 
				NodeList[np].bcode[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
			}
		}
	}
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
		if (!LoadCases[lcase].Read(Input, lcase, MODEX))
			return false;

	return true;
}


// Read element data
bool CDomain::ReadElements()
{
    EleGrpList = new CElementGroup[NUMEG];

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
        if (!EleGrpList[EleGrp].Read(Input))
            return false;
    
    return true;
}

//	Read history output message
bool CDomain::ReadHisMessage()
{
	
	Input >> Num_His_Output;
	His_freedom = new int[Num_His_Output * 2];
	for (int i = 0; i < Num_His_Output * 2; i++) {
		Input >> His_freedom[i];
	}

	return true;
}


//	Calculate column heights
void CDomain::CalculateColumnHeights()
{
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
        {
            CElement& Element = ElementGrp[Ele];

            // Generate location matrix
            Element.GenerateLocationMatrix();
            
            StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
			C_MassMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
        }
    }
    
    StiffnessMatrix->CalculateMaximumHalfBandwidth();
	C_MassMatrix->CalculateMaximumHalfBandwidth();
    
#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	Output->PrintColumnHeights();
#endif

}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();  //ElementGrp[0]->The first CElement object in ElementGroup(用ElementGroup类的对象控制CElement类的对象)
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
            Element.ElementStiffness(Matrix);
            StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
        }

		delete[] Matrix;
		Matrix = nullptr;
	}

#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	Output->PrintStiffnessMatrix();
#endif

}


//  Assemble the global mass matrix

void CDomain::AssembleMassMatrix()
{
	//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		CElementGroup& ElementGrp = EleGrpList[EleGrp];
		unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix(); // The length of consistent mass matrix is the same as the stiffness matrix
		double* Matrix = new double[size];

		//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
		{
			CElement& Element = ElementGrp[Ele];
			Element.ElementMass(Matrix);
			C_MassMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
		}
		delete[] Matrix;
		Matrix = nullptr;
	}

	COutputter* Outputter = COutputter::Instance();



#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	Output->PrintMassMatrix();
#endif
}

//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

//	Loop over for all concentrated loads in load case LoadCase
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        
        if(dof) // The DOF is activated
            Force[dof - 1] += LoadData->load[lnum];
	}

	return true;
}



//	Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//	and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
//	Allocate for global force/displacement vector
	Force = new double[NEQ];
    clear(Force, NEQ);

//  Create the banded stiffness matrix
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);

// Create the consistent banded mass matrix
	C_MassMatrix = new CSkylineMatrix<double>(NEQ);

// Create the lumped banded mass matrix
	L_MassMatrix = new double[NEQ];
	clear(L_MassMatrix, NEQ);

//	Calculate column heights
	CalculateColumnHeights();

//	Calculate address of diagonal elements in banded matrix
	StiffnessMatrix->CalculateDiagnoalAddress();
	C_MassMatrix->CalculateDiagnoalAddress();

//	Allocate for banded global stiffness matrix
    StiffnessMatrix->Allocate();
	C_MassMatrix->Allocate();

	COutputter* Output = COutputter::Instance();
	Output->OutputTotalSystemData();
}

