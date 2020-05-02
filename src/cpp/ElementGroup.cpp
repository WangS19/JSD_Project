/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "ElementGroup.h"
#include "Domain.h"
#include "Element.h"//yjl
//yjl
CElement::CElement() 
{
	NEN_ = 0;
	nodes_ = nullptr;
	ElementMaterial_ = nullptr;
	det_J = 0;
	XG[0][0] = 0;   XG[0][1] = -0.5773502691896;   XG[0][2] = -0.7745966692415;  XG[0][3] = -0.8611363115941;
	XG[1][0] = 0; 	XG[1][1] = 0.5773502691896;    XG[1][2] = 0;  				 XG[1][3] = -.3399810435849;
	XG[2][0] = 0; 	XG[2][1] = 0; 				   XG[2][2] = 0.7745966692415;   XG[2][3] = 0.3399810435849;
	XG[3][0] = 0; 	XG[3][1] = 0; 				   XG[3][2] = 0;  				 XG[3][3] = 0.8611363115941;
	WGT[0][0] = 2;  WGT[0][1] = 1;                 WGT[0][2] = 0.5555555555556;  WGT[0][3] = 0.3478548451375;
	WGT[1][0] = 0; 	WGT[1][1] = 1.0; 			   WGT[1][2] = 0.8888888888889;  WGT[1][3] = 0.6521451548625;
	WGT[2][0] = 0; 	WGT[2][1] = 0; 				   WGT[2][2] = 0.5555555555556;  WGT[2][3] = 0.6521451548625;
	WGT[3][0] = 0; 	WGT[3][1] = 0; 				   WGT[3][2] = 0;  				 WGT[3][3] = 0.3478548451375;
}




CNode* CElementGroup::NodeList_ = nullptr;

//! Constructor
CElementGroup::CElementGroup()
{
    if (!NodeList_)
    {
        CDomain* FEMData = CDomain::Instance();
        NodeList_ = FEMData->GetNodeList();
    }
    
    ElementType_ = ElementTypes::UNDEFINED;
    
    NUME_ = 0;
    ElementList_ = nullptr;
    
    NUMMAT_ = 0;
    MaterialList_ = nullptr;
}

//! Deconstructor
CElementGroup::~CElementGroup()
{
    if (ElementList_)
        delete [] ElementList_;
    
    if (MaterialList_)
        delete [] MaterialList_;
}

//! operator []
//! For the sake of efficiency, the index bounds are not checked
CElement& CElementGroup::operator[](unsigned int i)
{
    return *(CElement*)((std::size_t)(ElementList_) + i*ElementSize_);
}

//! Return index-th material in this element group
CMaterial& CElementGroup::GetMaterial(unsigned int index)
{
    return *(CMaterial*)((std::size_t)(MaterialList_) + index*MaterialSize_);
}

//! Calculate the size of the derived element and material class
void CElementGroup::CalculateMemberSize()
{
    switch (ElementType_)
    {
        case ElementTypes::UNDEFINED:
            std::cerr << "Setting element type to UNDEFINED." << std::endl;
            exit(5);
        case ElementTypes::Bar:
            ElementSize_ = sizeof(CBar);
            MaterialSize_ = sizeof(CBarMaterial);
            break;
		case ElementTypes::Q4:
			ElementSize_ = sizeof(CQ4);
			MaterialSize_ = sizeof(CQ4Material);
			break;
		case ElementTypes::AX8R:
			ElementSize_ = sizeof(CAX8R);
			MaterialSize_ = sizeof(CAX8RMaterial);
			break;//yjl
        default:
            std::cerr << "Type " << ElementType_ << " not available. See CElementGroup::CalculateMemberSize." << std::endl;
            exit(5);
            break;
    }
}

//! Allocate array of derived elements
void CElementGroup::AllocateElements(std::size_t size)
{
    switch(ElementType_)
    {
        case ElementTypes::Bar:
            ElementList_ = new CBar[size];
            break;
		case ElementTypes::Q4:
			ElementList_ = new CQ4[size];
			break;
		case ElementTypes::AX8R:
			ElementList_ = new CAX8R[size];
			break;
        default:
            std::cerr << "Type " << ElementType_ << " not available. See CElementGroup::AllocateElement." << std::endl;
            exit(5);
    }
}

//! Allocate array of derived materials
void CElementGroup::AllocateMaterials(std::size_t size)
{
    switch(ElementType_)
    {
        case ElementTypes::Bar:
            MaterialList_ = new CBarMaterial[size];
            break;
		case ElementTypes::Q4:
			MaterialList_ = new CQ4Material[size];
			break;
		case ElementTypes::AX8R://yjl
			MaterialList_ = new CAX8RMaterial[size];
			break;
        default:
            std::cerr << "Type " << ElementType_ << " not available. See CElementGroup::AllocateMaterial." << std::endl;
            exit(5);
    }
}

//! Read element group data from stream Input
bool CElementGroup::Read(ifstream& Input)
{
    Input >> (int&)ElementType_ >> NUME_ >> NUMMAT_;

	if (ElementType_ == 2) {
		Input >> N_G;   //Input the number of Gauss point in one direction
	}
	//yjl for CAX8R
	if (ElementType_ == 8) {
		Input >> N_G;   //Input the number of Gauss point in one direction
	}
    
    CalculateMemberSize();

    if (!ReadElementData(Input))
        return false;

    return true;
}

//  Read bar element data from the input data file
bool CElementGroup::ReadElementData(ifstream& Input)
{
//  Read material/section property lines
    AllocateMaterials(NUMMAT_);
    
//  Loop over for all material property sets in this element group
    for (unsigned int mset = 0; mset < NUMMAT_; mset++)
        if (!GetMaterial(mset).Read(Input, mset))
            return false;
    
//  Read element data lines
    AllocateElements(NUME_);

//  Get the number of Gauss point in one direction
	ElementList_->GetNG(N_G);
    
//  Loop over for all elements in this element group
    for (unsigned int Ele = 0; Ele < NUME_; Ele++)
        if (!(*this)[Ele].Read(Input, Ele, MaterialList_, NodeList_))
            return false;
    
    return true;
}
