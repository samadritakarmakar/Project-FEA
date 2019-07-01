#ifndef DIRICHLETBC_HPP
#define DIRICHLETBC_HPP
#include "TrialFunction.hpp"
#include "libGmshReader.h"
#include "Expression.hpp"
class DirichletBC
{
public:
    DirichletBC(TrialFunction &u, int PhysicalGroupNumber, umat& boolDiricletNodes):
        PhysclGrpNum (PhysicalGroupNumber), Msh(u.Msh),
        vctrLvlInternal(u.originalVctrLvl), currentNodeNumber(0)
    {
        DirichletBCNodes=std::vector<umat>(Msh->NumOfElementTypes);
        DirichletBCNodesFill=std::vector<umat>(Msh->NumOfElementTypes);
        NoOfNodes=std::vector<int>(Msh->NumOfElementTypes);
        NodePositions=std::vector<umat>(Msh->NumOfElementTypes);
        ExprssnPositions=std::vector<umat>(Msh->NumOfElementTypes);
        for (int ElementType=0; ElementType<Msh->NumOfElementTypes; ElementType++)
        {
            DirichletBCNodes[ElementType]=vectorise(unique(Msh->ElmntPhysclGrpNodes[ElementType][PhysclGrpNum]));
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
     void x(vec &X)
    {
         mat Coordinates;
         for (int ElementType=0; ElementType<Msh->NumOfElementTypes; ElementType++)
         {
             Coordinates=Msh->NodalCoordinates.rows(DirichletBCNodes[ElementType]);
         }
        // X= vectorise(Coordinates(currentNodeNumber,span(0,vctrLvlInternal-1)));
         X= vectorise(Coordinates.cols(0,vctrLvlInternal-1));
     }

     /// Set the Expression used to evaluate the value of Dirichlet Boundary at each node
     void SetDirichletBCExpression(Expression& DeclaredExprssnClss)
         {
             ExprssnInternal= &DeclaredExprssnClss;
         }

     /// Applies the the given Dirichlet BC.
     void ApplyBC(sp_mat &A, mat& b)
     {
         cout<<"Applying Dirichlet BC over Physical Group "<<Msh->PhysicalGroupName[PhysclGrpNum]<<"\n";
         // Below the The rows and cols of matrix 'A' with dirichlet BC are set to zero.
         // The diagonal terms with  dirichlet BC are set to 1.
         // The rows of vector 'b' are set the with Diriclet values.
         for (int ElementType=0; ElementType<Msh->NumOfElementTypes; ElementType++)
         {
             umat locations_1s;
             umat uniquePositions;
             uniquePositions.set_size(1,0);
             vec XTemp;
             x(XTemp);
             for (int NodeNumber=0; NodeNumber<NoOfNodes[ElementType]; NodeNumber++)
             {
                 currentNodeNumber=NodeNumber;
                 umat Positions1=DirichletBCNodesFill[ElementType].row(NodeNumber);
                 uniquePositions=join_horiz(uniquePositions, Positions1);
                 //vec X;
                 //x(X);
                 vec X=XTemp.row(NodeNumber);
                 vec Expressn= ExprssnInternal->Eval(X);
                 b.rows(Positions1)=Expressn.rows(ExprssnPositions[ElementType]);
                 //cout<<"DirchletBC Applied at, "<<Msh->NodalCoordinates.row(Positions1(0,0)/vctrLvlInternal)<<"\n";
                 /*sp_mat Azero=speye(A.n_rows, A.n_cols);
                 for (int MatPosition=0;MatPosition<Positions1.n_cols;MatPosition++)
                 {
                     //A.row(Positions1(0,MatPosition))=zeros(1,A.n_cols);
                     //A.col(Positions1(0,MatPosition))=zeros(A.n_rows,1);
                     Azero(Positions1(0,MatPosition),Positions1(0,MatPosition))=0.0;
                     Azero=sp_mat(Azero);
                 }
                 A=A*Azero;
                 A=(A.t()*Azero).t();*/
             }
             locations_1s.set_size(2, uniquePositions.n_cols);
             locations_1s=repmat(uniquePositions,2,1);
             //This saves a set of sparse matrix whose all elements are zero except where Dirchlet BC conditons are applied.
             //Where Dirichlet BC Conditions are applied, diagonal elements are '1.0';
             sp_mat ATemp=sp_mat(false, locations_1s, ones(locations_1s.n_cols,1), A.n_rows, A.n_cols, true, false);
             //---------------------------
             sp_mat Azero=speye(A.n_rows, A.n_cols);
             Azero-=ATemp;//Where Dirichlet BC Conditions are applied, diagonal elements are '0.0' other diagonals are 1.0;
             A=A*Azero; //Sets Columns corresponding to Dirichlet BC Conditions to zero.
             A=(A.t()*Azero).t(); //Sets Rows corresponding to Dirichlet BC Conditions to zero.
             //--------------------------
             //cout<<mat(ATemp)<<"\n";
             A=A+ATemp; //Where Dirichlet BC Conditions are applied, diagonal elements are '1.0'.
                        //Corresponding Rows and columns are zero;
             cout<<"Done Applying Dirichlet BC!!!\n";
         }
     }
private:
    //TrialFunction& u_Internal;
    int& PhysclGrpNum;
    int& vctrLvlInternal;
    int currentNodeNumber;
    libGmshReader::MeshReader *Msh;
    std::vector<umat> DirichletBCNodes, DirichletBCNodesFill;
    std::vector<int> NoOfNodes;
    std::vector<umat> NodePositions;
    std::vector<umat> ExprssnPositions;
    Expression* ExprssnInternal;
    mat Coordinates;
};

#endif // DIRICHLETBC_HPP
