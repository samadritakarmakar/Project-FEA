#ifndef DIRICHLETBC_HPP
#define DIRICHLETBC_HPP
#include "TrialFunction.hpp"
#include "libGmshReader.h"
#include "Expression.hpp"
class DirichletBC
{
public:
    DirichletBC(TrialFunction &u, int PhysicalGroupNumber, umat& boolDiricletNodes):
        u_Internal(u), PhysclGrpNum (PhysicalGroupNumber), Msh(u.Msh),
        vctrLvlInternal(u.originalVctrLvl), currentNodeNumber(0)
    {
        //ExprssnInternal = std::shared_ptr<Expression>(new Expression(u_Internal.originalVctrLvl));
        DirichletBCNodes=std::vector<umat>(Msh->NumOfElementTypes);
        DirichletBCNodesFill=std::vector<umat>(Msh->NumOfElementTypes);
        NoOfNodes=std::vector<int>(Msh->NumOfElementTypes);
        NodePositions=std::vector<umat>(Msh->NumOfElementTypes);
        ExprssnPositions=std::vector<umat>(Msh->NumOfElementTypes);
        for (int ElementType=0; ElementType<Msh->NumOfElementTypes; ElementType++)
        {
            DirichletBCNodes[ElementType]=unique(Msh->ElmntPhysclGrpNodes[ElementType][PhysclGrpNum]);
            NoOfNodes[ElementType]=DirichletBCNodes[ElementType].n_rows;
            GetNodePostions(NodePositions[ElementType],DirichletBCNodes[ElementType],vctrLvlInternal);
            DirichletBCNodesFill[ElementType].set_size(NodePositions[ElementType].n_rows,1);
            int DiricletNodeSize=1;
            int boolDiriclet=0;
            for (int colNo=0;colNo<NodePositions[ElementType].n_cols;colNo++)
            {
                if(bool(boolDiricletNodes(0,boolDiriclet)))
                {
                    DirichletBCNodesFill[ElementType].resize(NodePositions[ElementType].n_rows, DiricletNodeSize);
                    DirichletBCNodesFill[ElementType].col(DiricletNodeSize-1)=NodePositions[ElementType].col(colNo);
                    DiricletNodeSize++;
                }
                boolDiriclet++;
                if(boolDiriclet>=vctrLvlInternal)
                {
                    boolDiriclet=0;
                }
            }
            boolDiriclet=0;
            int ExprssnSize=1;
            ExprssnPositions[ElementType].set_size(1, ExprssnSize);
            for (int ExprssIndx=0; ExprssIndx<vctrLvlInternal; ExprssIndx++)
            {
                if(bool(boolDiricletNodes(0,ExprssIndx)))
                {
                    //cout<<"ExprssIndx = "<<ExprssIndx<<"\n";
                    ExprssnPositions[ElementType].resize(1,ExprssnSize);
                    ExprssnPositions[ElementType](0,ExprssnSize-1)=ExprssIndx;
                    ExprssnSize++;
                }
                boolDiriclet++;
            }
        }
        //cout<<NodePositions[0];
        //cout<<DirichletBCNodesFill[0];
        //cout<<ExprssnPositions[0];
        //cout<<"Diriclet Nodes =\n"<<DirichletBCNodes[0];
    }

    /*virtual mat Expression(mat& x)
    {
        mat Exprssn=zeros(vctrLvlInternal,1);
        return Exprssn;
    }*/

    ///Returns a matrix of coordinates [x, y, z] of all current node number
    /// within the Physical Enitity.
     vec x()
    {
         mat Coordinates;
         for (int ElementType=0; ElementType<Msh->NumOfElementTypes; ElementType++)
         {
             Coordinates=Msh->NodalCoordinates.rows(DirichletBCNodes[ElementType]);
         }
         return vectorise(Coordinates(currentNodeNumber,span(0,vctrLvlInternal-1)));
     }

     void SetDirichletBC(Expression& DeclaredExprssnClss)
         {
             ExprssnInternal= &DeclaredExprssnClss;
         }

     void ApplyBC(sp_mat &A, mat& b)
     {
         for (int ElementType=0; ElementType<Msh->NumOfElementTypes; ElementType++)
         {
             //umat NodePosition;
             //GetNodePostions(NodePosition, DirichletBCNodes[ElementType], vctrLvlInternal);
             for (int NodeNumber=0; NodeNumber<NoOfNodes[ElementType]; NodeNumber++)
             {
                 currentNodeNumber=NodeNumber;
                 //umat Positions1=NodePosition.row(NodeNumber);
                 umat Positions1=DirichletBCNodesFill[ElementType].row(NodeNumber);
                 //mat Expressn=Expression(x);
                 vec X=x();
                 vec Expressn= ExprssnInternal->Eval(X);
                 b.rows(Positions1)=Expressn.rows(ExprssnPositions[ElementType]);
                 for (int MatPosition=0;MatPosition<Positions1.n_cols;MatPosition++)
                 {
                     A.row(Positions1(0,MatPosition))=zeros(1,A.n_cols);
                     A.col(Positions1(0,MatPosition))=zeros(A.n_rows,1);
                     A(Positions1(0,MatPosition),Positions1(0,MatPosition))=1;
                 }
             }
         }
     }
private:
    TrialFunction& u_Internal;
    int& PhysclGrpNum;
    int& vctrLvlInternal;
    int currentNodeNumber;
    libGmshReader::MeshReader *Msh;
    std::vector<umat> DirichletBCNodes, DirichletBCNodesFill;
    std::vector<int> NoOfNodes;
    std::vector<umat> NodePositions;
    std::vector<umat> ExprssnPositions;
    //std::shared_ptr<Expression> ExprssnInternal;
    Expression* ExprssnInternal;
};

#endif // DIRICHLETBC_HPP
