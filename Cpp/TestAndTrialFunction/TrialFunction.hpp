#ifndef TRIALFUNCTION_HPP
#define TRIALFUNCTION_HPP
#include "libGmshReader.h"
#include "FEMtools.h"
#include "LagrangeShapeFunction.hpp"
#include "Jacobian.hpp"
#include <vector>
#include <armadillo>
#include <assert.h>
using namespace arma;
class OtherData
{
public:
    int ElementType=0;
    //int GuassPntCounter=0;
    //libGmshReader::MeshReader Mesh;
};
class TrialFunction
{
public:
    mat (TrialFunction::*Pntr_Calc_N)(mat , OtherData);
    std::vector<int> NoOfGaussPts;
    std::vector<int> NoOfElements;
    std::vector<umat> ElmntNodes;
    libGmshReader::MeshReader *Msh;
    int MeshDimension;
    int vectorLvl;
    int originalVctrLvl;
    int NoOfElementTypes;

    /// To be used only for declaring, must not forget to Initialize.
    TrialFunction()
    {
        //---Do--Nothing
    }
    TrialFunction(libGmshReader::MeshReader& Mesh, int MeshDimension, int vectorLevel)
    {
        Msh =new libGmshReader::MeshReader(Mesh, MeshDimension);
        NewMeshInstanceCreated=true;
        //std::cout<<"Constructor runs file name type !!\n";
        SetDimension();
        SetCommonVariables(vectorLevel);
        if(vectorLevel<MeshDimension && vectorLevel!=1)
        {
            vectorLvl=1;
        }
        SetDimension();
        GetNumberOfVariables();
        SetGaussPtBasedVariables();
        // Keep derivatives of Shape Functions w.rt. points within the reference elements ready.
        Generate_dN_by_dEps ();
    }

    /// Must use this to initialize the Trial Function
    TrialFunction(libGmshReader::MeshReader& Mesh, int& vectorLevel)
    {
        NewMeshInstanceCreated=false;
        //std::cout<<"Constructor runs Mesh type !!\n";
        SetMesh(Mesh);
        SetCommonVariables(vectorLevel);
        SetDimension();
        CheckVectorLevel();
        GetNumberOfVariables();
        SetGaussPtBasedVariables();
        // Keep derivatives of Shape Functions w.rt. points within the reference elements ready.
        Generate_dN_by_dEps ();
        /*std::vector<libGmshReader::MeshReader> MeshLowerDim (MeshDimension-1);

        int dim;
        int vectorLvlNew=vectorLvl;
        for (dim=MeshDimension-1; dim>=2; dim--)
        {
            if(vectorLvlNew!=1)
            {
                vectorLvlNew--;
            }
            cout<<"vectorLvlNew = "<<vectorLvlNew<<" dim= "<<dim<<"\n";
            MeshLowerDim[dim-1]=libGmshReader::MeshReader(Msh->ElementData::fileName, dim);
            //static TrialFunction u_down(MeshLowerDim[dim],vectorLvlNew, dim);
        }*/
    }

    /*TrialFunction(libGmshReader::MeshReader& Mesh, int& vectorLevel, int& Dimension)
    {
        NewMeshInstanceCreated=false;
        SetMesh(Mesh);
        SetCommonVariables(vectorLevel);
        SetDimension(Dimension);
        CheckVectorLevel();
        GetNumberOfVariables();
        SetGaussPtBasedVariables();
        // Keep derivatives of Shape Functions w.rt. points within the reference elements ready.
        Generate_dN_by_dEps ();
    }
*/
    void SetMesh(libGmshReader::MeshReader& Mesh)
    {
        Msh=&Mesh;
    }


    void SetCommonVariables(int& vectorLevel)
    {
        vectorLvl=vectorLevel;
        originalVctrLvl=vectorLevel;
        N=std::vector<mat> (Msh->NumOfElementTypes);
        Phi=std::vector<LagrangeShapeFunction> (Msh->NumOfElementTypes);
        GaussData=std::vector<mat> (Msh->NumOfElementTypes);
        Weight=std::vector<mat> (Msh->NumOfElementTypes);
        GaussPointx=std::vector<mat> (Msh->NumOfElementTypes);
        GaussPointy=std::vector<mat> (Msh->NumOfElementTypes);
        GaussPointz=std::vector<mat> (Msh->NumOfElementTypes);
        NoOfGaussPts=std::vector<int> (Msh->NumOfElementTypes);
        u=std::vector<std::vector<sp_mat>>(Msh->NumOfElementTypes);
        dN_by_dEps=std::vector<std::vector<mat>>(Msh->NumOfElementTypes);
        NoOfElements=std::vector<int> (Msh->NumOfElementTypes);
        ElmntNodes=std::vector<umat> (Msh->NumOfElementTypes);
        NoOfElementTypes=Msh->NumOfElementTypes;
    }

    void SetDimension()
    {
        MeshDimension=Msh->ElementData::dim;
    }

    void SetDimension(int& Dimension)
    {
        MeshDimension=Dimension;
    }

    void CheckVectorLevel()
    {
        if (vectorLvl!=1 && vectorLvl!=MeshDimension)
        {
            std::cout<<"Dimension ("<<MeshDimension<<") must match the Vector Level "<<vectorLvl<<"\n";
            std::cout<<"or Vector Level should be equal to 1!!!\n";
            throw;
        }
    }

    void SetGaussPtBasedVariables()
    {
        for (int ElementType = 0; ElementType<Msh->NumOfElementTypes; ++ElementType)
        {
            Generate_GaussPoints_Weights_ShapeFunctions(ElementType);
            NoOfGaussPts[ElementType]=GaussPointx[ElementType].n_rows;
            u[ElementType]=std::vector<sp_mat>(NoOfGaussPts[ElementType]);
            dN_by_dEps[ElementType]=std::vector<mat>(NoOfGaussPts[ElementType]);
            /*N[ElementType]=Phi[ElementType].GetShapeFunction(GaussPointx[ElementType],GaussPointy[ElementType],
                                                             GaussPointz[ElementType]);*/
            for (int GaussPt=0;GaussPt<NoOfGaussPts[ElementType];GaussPt++)
            {
                mat Ncol=N[ElementType].col(GaussPt);
                u[ElementType][GaussPt]=vectorizeQuantity(Ncol,vectorLvl);
                //cout<<mat(u[ElementType][GaussPt])<<"\n";
            }
        }
    }

    void Regenerate_u()
    {
        for (int ElementType = 0; ElementType<Msh->NumOfElementTypes; ++ElementType)
        {
            u[ElementType]=std::vector<sp_mat>(NoOfGaussPts[ElementType]);
            for (int GaussPt=0;GaussPt<NoOfGaussPts[ElementType];GaussPt++)
            {
                mat Ncol=N[ElementType].col(GaussPt);
                u[ElementType][GaussPt]=vectorizeQuantity(Ncol,originalVctrLvl);
                //cout<<mat(u[ElementType][GaussPt])<<"\n";
            }
        }
    }

    /// Sets the NoOfElements and ElmntNodes or 'Connectivity Matrix' for a certain Element Type.
    void GetNumberOfVariables()
    {
        for (int ElementType = 0; ElementType<Msh->NumOfElementTypes; ++ElementType)
        {
            NoOfElements[ElementType]= Msh->ElementNodes[ElementType].n_rows;
            ElmntNodes[ElementType]= Msh->ElementNodes[ElementType];
        }
    }


    /// Destructor
    ~TrialFunction()
    {
        //std::cout<<"Destructor runs!!\n";
        //delete [] u_down;
        if (NewMeshInstanceCreated)
        {
            //std::cout<<"Mesh Destroyed of dim"<<Msh->ElementData::dim <<"!\n";
            //delete []Msh;
        }
    }

    /* sp_mat dot_vectrLvl_grad_u(vec a, sp_mat grad_u)
    {
        if (vectorLvl==1)
        {
            //mat aMatrx=repmat(a.t(),grad_u.n_rows,1);
            //return sp_mat(aMatrx%grad_u);
            mat aMatrx=a.rows(1,MeshDimension).t();
            //cout<<grad_u;
            return sp_mat(aMatrx*grad_u);
        }
        else
        {
            mat vctr=a(span(0,vectorLvl-1)).t();
            sp_mat vecMatrx(vctr.n_cols,vectorLvl*vctr.n_cols);
            for (int i=0; i<vectorLvl; i++)
            {
                //cout<<"col ="<<i*vectorLvl<<":"<<i*vectorLvl+vectorLvl-1;
                vecMatrx.cols(i*vectorLvl,i*vectorLvl+vectorLvl-1)=sp_mat(diagmat(vctr.t()));
            }
            return vecMatrx*grad_u;
        }
    }
    */
    /// Gets the interpolated value of [x,y,z] at a cetrian Gauss point
    mat Get_x(int ElementType, int ElementNumber, int GaussPntr)
    {
        umat NodesAtElmntNmbr=Msh->ElementNodes[ElementType].row(ElementNumber);
        mat Coordinates=Msh->NodalCoordinates.rows(NodesAtElmntNmbr);
        mat x=Coordinates.cols(0,2).t()*N[ElementType].col(GaussPntr);
        return x;
    }



    /// Gives an output for trace of grad u
     mat Get_trace_grad_u(int ElementType, int ElementNumber, int GaussPntr)
    {
        mat dN_by_dx_atGaussPt=Get_dN_by_dx(ElementType, ElementNumber, GaussPntr);
        int cols_dN_dx=Get_dN_by_dx(ElementType,0,0).n_cols;
        mat trace=repmat(vectorise(dN_by_dx_atGaussPt,1),vectorLvl*cols_dN_dx,1);
        return trace;
    }

    /// Returns Kronecker Delta that can be used with gradients on u.
     sp_mat Get_KroneckerDelta(int ElementType)
    {
        int GaussPntr=0;
        int rows_dN_dx=Get_dN_by_dx(ElementType,0,0).n_rows;
        int cols_dN_dx=Get_dN_by_dx(ElementType,0,0).n_cols;
        sp_mat KroneckerDelta(vectorLvl*cols_dN_dx,vectorLvl*rows_dN_dx);
        for (int row=0; row<cols_dN_dx; row++)
        {
            for (int col=0; col<rows_dN_dx; col++)
            {
                KroneckerDelta(row*cols_dN_dx+row,col*vectorLvl+row)=1;
            }
        }
        return KroneckerDelta;
    }

    /// Returns curl of u.
     sp_mat Get_curl_u(int ElementType, int ElementNumber, int GaussPntr)
    {
        mat dN_by_dx_atGaussPt=Get_dN_by_dx(ElementType, ElementNumber, GaussPntr);
        int rows_dN_dx=dN_by_dx_atGaussPt.n_rows;
        int cols_dN_dx=dN_by_dx_atGaussPt.n_cols;
        sp_mat curl_u(vectorLvl, vectorLvl*cols_dN_dx);
        if (vectorLvl==3)
        {
            for (int row = 0; row < cols_dN_dx; ++row)
            {
                int colAddStart=vectorLvl*row;
                int colAddEnd=vectorLvl*row+(vectorLvl-1);
                curl_u(2,colAddStart+1)=dN_by_dx_atGaussPt(row,0);
                curl_u(0,colAddStart+2)=dN_by_dx_atGaussPt(row,1);
                curl_u(1,colAddStart+0)=dN_by_dx_atGaussPt(row,2);
                curl_u(1,colAddStart+2)=-dN_by_dx_atGaussPt(row,0);
                curl_u(2,colAddStart+0)=-dN_by_dx_atGaussPt(row,1);
                curl_u(0,colAddStart+1)=-dN_by_dx_atGaussPt(row,2);
            }
        }
        else
        {
            std::cout<<"Curl of Only VectorLevel 3 is supported!!!\n";
            throw;
        }
        return curl_u;
    }

    /// Returns the Gradient of u for a particular gauss point of a particular
    /// element which is a particular ElementType
     sp_mat Get_grad_u(int ElementType, int ElementNumber, int GaussPntr)
    {
        mat dN_by_dx_atGaussPt=Get_dN_by_dx(ElementType, ElementNumber, GaussPntr);
        int rows_dN_dx=dN_by_dx_atGaussPt.n_rows;
        int cols_dN_dx=dN_by_dx_atGaussPt.n_cols;
        sp_mat grad_u(vectorLvl*rows_dN_dx, vectorLvl*cols_dN_dx);
        int row=0;
        //cout<<"vectorLvl ="<<vectorLvl<<"\n";
        for (int cntr=0; cntr<rows_dN_dx; cntr++)
        {
            for (int col=0; col<vectorLvl; col++)
            {
                grad_u(row,span(col*vectorLvl, col*vectorLvl+cols_dN_dx-1))=dN_by_dx_atGaussPt.row(cntr);
                row++;
            }
        }
        return grad_u.t();
    }
    /// Returns the Transpose of Gradient of u for a particular gauss point of a particular
    /// element which is a particular ElementType
     sp_mat GetTranspose_grad_u(int ElementType, int ElementNumber, int GaussPntr)
    {
        mat dN_by_dx_atGaussPt=Get_dN_by_dx(ElementType, ElementNumber, GaussPntr);
        //cout<<"dN_by_dx_atGaussPt=\n"<<dN_by_dx_atGaussPt<<"\n";
        int rows_dN_dx=dN_by_dx_atGaussPt.n_rows;
        int cols_dN_dx=dN_by_dx_atGaussPt.n_cols;
        sp_mat Transpose_grad_u(vectorLvl*cols_dN_dx,vectorLvl*rows_dN_dx);
        for (int row=0; row<rows_dN_dx; row++)
        {
            for (int col=0; col<cols_dN_dx; col++)
            {
                //cout<<"row= "<<col*vectorLvl<<","<<col*vectorLvl+vectorLvl-1<<
                //      "; col= "<<row*vectorLvl<<","<<row*vectorLvl+vectorLvl-1<<"\n";
                Transpose_grad_u(span(col*vectorLvl,col*vectorLvl+vectorLvl-1),span(row*vectorLvl,row*vectorLvl+vectorLvl-1))=
                        dN_by_dx_atGaussPt(row,col)*speye(vectorLvl,vectorLvl);
            }
        }
        return Transpose_grad_u;
    }

    /// Returns the Gradient of N (dN_by_dx) for a particular gauss point of a particular
    /// element which is a particular ElementType
     mat Get_dN_by_dx(int ElementType, int ElementNumber, int GaussPntr)
    {
        mat F;
        Get_F(ElementType, ElementNumber, GaussPntr, F);
        mat dN_by_dx=dN_by_dEps[ElementType][GaussPntr]*inv(F);
        return dN_by_dx;
    }

    /// This function generates the F matrix. This in mathematical terms is the Jacobian dx/dEps
    void Get_F(int ElementType, int ElementNumber, int GaussPntr, mat& F)
    {
        umat NodesAtElmntNmbr=Msh->ElementNodes[ElementType].row(ElementNumber);
        //cout<<"Gmsh Node Tags are "<<Msh->GmshNodeTag[ElementType].row(ElementNumber);
        //cout<<"ElementNodes =\n"<<NodesAtElmntNmbr;
        //cout<<"Coodinates of "<<ElementNumber<<" are \n"<<Msh->NodalCoordinates.rows(NodesAtElmntNmbr);
        mat Coordinates=Msh->NodalCoordinates.rows(NodesAtElmntNmbr);
        //cout<<"Coodinates of "<<ElementNumber<<" are \n"<<Coordinates.cols(0,MeshDimension-1)<<"\n";
        //coords of x for dim 1; x & y for dim 2; x, y & z for dim 3;
        F=Coordinates.cols(0,MeshDimension-1).t()*dN_by_dEps[ElementType][GaussPntr];
        //cout<<"Coordinates.cols(0,MeshDimension-1).t()=\n"<<Coordinates.cols(0,MeshDimension-1).t();
        //cout<<"\ndN_by_dEps[ElementType][GaussPntr];\n"<<dN_by_dEps[ElementType][GaussPntr];
    }

    /// This function generates the Jacobian of dN_by_dEps
    void Generate_dN_by_dEps ()
    {
        OtherData Data;
        int dim=MeshDimension;
        mat x;
        x.set_size(1,dim);
        for (int ElmntTypCountr=0; ElmntTypCountr<Msh->NumOfElementTypes; ElmntTypCountr++)
        {
            int &i=ElmntTypCountr;
            Data.ElementType=ElmntTypCountr;
            for (int GaussPntr=0; GaussPntr<GaussPointx[i].n_rows; GaussPntr++)
            {
                if(MeshDimension==1)
                {
                    x=GaussPointx[i].row(GaussPntr);
                }
                else if (MeshDimension==2)
                {
                    x=join_vert(GaussPointx[i].row(GaussPntr), GaussPointy[i].row(GaussPntr));
                }
                else if (MeshDimension==3)
                {
                    mat x1=join_vert(GaussPointx[i].row(GaussPntr), GaussPointy[i].row(GaussPntr));
                    x=join_vert(x1,GaussPointz[i].row(GaussPntr));
                }
                Pntr_Calc_N=&TrialFunction::Calculate_N_GaussPointWise;
                //Data.GuassPntCounter=GaussPntr;
                dN_by_dEps[i][GaussPntr]= Jacobian(this,(this->Pntr_Calc_N), x, Data);
                //cout<<"dN_by_dEps at Gauss Pnt "<<GaussPntr<<" =\n"<<dN_by_dEps[i][GaussPntr];
            }
        }
    }
    /// This function calculates the shape function for a gauss point of an element,
    /// belonging to a particular element type.
     mat Calculate_N_GaussPointWise(mat x1, OtherData Data)
    {
        vec x=vectorise(x1);
        //int &ElemntTyp=Data.ElementType;
        //int &GaussPt=Data.GuassPntCounter;
        mat GaussPtx,GaussPty,GaussPtz, N_AtGaussPt;
        if (MeshDimension==1)
        {
            //cout<<"x= "<<x<<"\n";
            GaussPtx=x.row(0);
            N_AtGaussPt=Phi[Data.ElementType].GetShapeFunction(GaussPtx);
        }
        else if (MeshDimension==2)
        {
            //cout<<"x= "<<x<<"\n";
            GaussPtx=x.row(0);
            GaussPty=x.row(1);
            N_AtGaussPt=Phi[Data.ElementType].GetShapeFunction(GaussPtx,GaussPty);
        }
        else if(MeshDimension==3)
        {
            //cout<<"x= "<<x<<"\n";
            GaussPtx=x.row(0);
            GaussPty=x.row(1);
            GaussPtz=x.row(2);
            N_AtGaussPt=Phi[Data.ElementType].GetShapeFunction(GaussPtx,GaussPty,GaussPtz);
        }
        //sp_mat u_AtGaussPt=vectorizeQuantity(N_AtGaussPt, vectorLvl);
        return N_AtGaussPt;
    }

     double WeightAt(int ElementType, int GaussPntr)
    {
        return Weight[ElementType](GaussPntr);
    }

    /// Returns the Gauss Weights for a particular ElementType
     mat GetWeights(int ElementType)
    {
        return Weight[ElementType];
    }
    /// Returns the Gauss Points of x-coordinates for a particular ElementType
     mat GetGaussPoints_x(int ElementType)
    {
        if (MeshDimension>=1)
            return GaussPointx[ElementType];
        else
        {
            cout<<"Dimension has to be at least 1 or Greater to return GaussPoints of x";
            throw;
            return {0};
        }
    }

    /// Returns the Gauss Points of y-coordinates for a particular ElementType
     mat GetGaussPoints_y(int ElementType)
    {
        if (MeshDimension>=2)
            return GaussPointy[ElementType];
        else
        {
            cout<<"Dimension has to be at least 2 or Greater to return GaussPoints of y";
            throw;
            return {0};
        }
    }
    /// Returns the Gauss Points of z-coordinates for a particular ElementType
     mat GetGaussPoints_z(int ElementType)
    {
        if (MeshDimension>=3)
            return GaussPointy[ElementType];
        else
        {
            cout<<"Dimension has to be at least 2 or Greater to return GaussPoints of z";
            throw;
            return {0};
        }
    }

    /// Returns the Number of Gauss Points per Element Type.
     int GetNumberOfGaussPoints(int ElementType)
    {
        return NoOfGaussPts[ElementType];
    }
    /// Returns the matrix of the Shape function of u.
     sp_mat Get_u(int ElementType, int GaussPntr)
    {
        return u[ElementType][GaussPntr];
    }
protected:
    std::vector<mat> GaussData;
    std::vector<mat> Weight;
    std::vector<mat> GaussPointx;
    std::vector<mat> GaussPointy;
    std::vector<mat> GaussPointz;
    std::vector<std::vector<sp_mat>> u;
    std::vector<std::vector<mat>> dN_by_dEps;
    std::vector<LagrangeShapeFunction> Phi;
    std::vector<mat> N;
    bool NewMeshInstanceCreated=false;
    //std::string MshFileName;

    //TrialFunction *u_down;

    /// 'Variable' must be in form of a vector
    /// This function builds the representation similar to
    /// how shape functions are represented generally in text books
    /// [N1 0 0; 0 N1 0; 0 0 N1; N2 0 0; 0 N2 0; 0 0 N2; N3 0 0; 0 N3 0; 0 0 N3]
     sp_mat vectorizeQuantity(mat& Variable, int VectorLevel)
    {
        vec VariableVector=vectorise(Variable);
        int NoOfRows=VariableVector.n_rows;
        sp_mat I=speye(VectorLevel, VectorLevel);
        //cout<<"I=\n"<<mat(I)<<"\n";
        //cout<<"N=\n"<<Variable<<"\n";
        sp_mat Matrx(VectorLevel, NoOfRows*VectorLevel);
        for (int i=0;i<NoOfRows;i++)
        {
            Matrx.cols(VectorLevel*i,VectorLevel*i+(VectorLevel-1))=VariableVector(i)*I;
        }
        return Matrx;
    }

    void CalculateGaussPointsAndWeights(mat &wt, mat& GaussPntx, mat& GaussPnty, int ElementType)
    {
        wt.set_size(pow(GaussData[ElementType].n_rows,2),1);
        //cout<<"Wt size =("<<wt.n_rows<<","<<wt.n_cols<<")\n";
        GaussPntx.set_size(wt.n_rows,1);
        GaussPnty.set_size(wt.n_rows,1);
        int pos=0;
        for(int i=0; i<GaussData[ElementType].n_rows; i++)
        {
            for (int j=0; j<GaussData[ElementType].n_rows; j++)
            {
                wt(pos,0)=GaussData[ElementType](i,1)*GaussData[ElementType](j,1);
                GaussPntx(pos,0)=GaussData[ElementType](i,2);
                GaussPnty(pos,0)=GaussData[ElementType](j,2);
                //cout<<"i= "<<i<<" j= "<<j<<" pos="<<pos<<"\n";
                pos++;
            }
        }
    }

    void CalculateGaussPointsAndWeights(mat &wt, mat& GaussPntx, mat& GaussPnty, mat& GaussPntz, int ElementType)
    {
        wt.set_size(pow(GaussData[ElementType].n_rows,3),1);
        GaussPntx.set_size(wt.n_rows,1);
        GaussPnty.set_size(wt.n_rows,1);
        GaussPntz.set_size(wt.n_rows,1);
        int pos=0;
        for(int i=0; i<GaussData[ElementType].n_rows; i++)
        {
            for (int j=0; j<GaussData[ElementType].n_rows; j++)
            {
                for (int k=0; k<GaussData[ElementType].n_rows; k++)
                {
                    wt(pos,0)=GaussData[ElementType](i,1)*GaussData[ElementType](j,1)*GaussData[ElementType](k,1);
                    GaussPntx(pos,0)=GaussData[ElementType](i,2);
                    GaussPnty(pos,0)=GaussData[ElementType](j,2);
                    GaussPntz(pos,0)=GaussData[ElementType](k,2);
                    pos++;
                }
            }
        }
    }

    void Generate_GaussPoints_Weights_ShapeFunctions(int ElementType)
    {
        Phi[ElementType]=LagrangeShapeFunction(*Msh,ElementType);
        GaussData[ElementType].load(FEMtools::LoadGaussFile(*Msh,ElementType));
        if(Msh->GmshElementNameOnly[ElementType].compare("Line")==0)
        {
            Weight[ElementType]=GaussData[ElementType].col(1);
            GaussPointx[ElementType]=GaussData[ElementType].col(2);
            N[ElementType]=Phi[ElementType].GetShapeFunction(GaussPointx[ElementType]);
        }
        else if(Msh->GmshElementNameOnly[ElementType].compare("Triangle")==0)
        {
            Weight[ElementType]=GaussData[ElementType].col(1);
            GaussPointx[ElementType]=GaussData[ElementType].col(2);
            GaussPointy[ElementType]=GaussData[ElementType].col(3);
            N[ElementType]=Phi[ElementType].GetShapeFunction(GaussPointx[ElementType],GaussPointy[ElementType]);
        }
        else if (Msh->GmshElementNameOnly[ElementType].compare("Quadrilateral")==0)
        {
            CalculateGaussPointsAndWeights(Weight[ElementType],GaussPointx[ElementType],GaussPointy[ElementType],ElementType);
            N[ElementType]=Phi[ElementType].GetShapeFunction(GaussPointx[ElementType],GaussPointy[ElementType]);
        }
        else if (Msh->GmshElementNameOnly[ElementType].compare("Tetrahedron")==0)
        {
            Weight[ElementType]=GaussData[ElementType].col(1);
            GaussPointx[ElementType]=GaussData[ElementType].col(2);
            GaussPointy[ElementType]=GaussData[ElementType].col(3);
            GaussPointz[ElementType]=GaussData[ElementType].col(4);
            N[ElementType]=Phi[ElementType].GetShapeFunction(GaussPointx[ElementType],GaussPointy[ElementType],
                                                             GaussPointz[ElementType]);
        }
        else if(Msh->GmshElementNameOnly[ElementType].compare("Hexahedron")==0
                ||Msh->GmshElementNameOnly[ElementType].compare("Prism")==0
                ||Msh->GmshElementNameOnly[ElementType].compare("Pyramid")==0 )
        {
            CalculateGaussPointsAndWeights(Weight[ElementType],GaussPointx[ElementType],GaussPointy[ElementType],
                                           GaussPointz[ElementType],ElementType);
            N[ElementType]=Phi[ElementType].GetShapeFunction(GaussPointx[ElementType],GaussPointy[ElementType],
                                                             GaussPointz[ElementType]);
        }
    }

    /// This function is strictly for debugging
    void PrintAll_dN_by_dEps()
    {
        for (int ElementType=0; ElementType<Msh->NumOfElementTypes; ElementType++)
        {
            for (int GaussPntr=0; GaussPntr<GaussPointx[ElementType].n_rows; GaussPntr++)
            {
                cout<<"dN_by_dEps at Element Type= "<<ElementType<<" GaussPntr= "<<GaussPntr<<"\n";
                cout<<"Where Gaussx ="<<GaussPointx[ElementType](GaussPntr,0)<< "Gauss y ="<<GaussPointy[ElementType](GaussPntr,0)<<"Gauss z ="<<GaussPointz[ElementType](GaussPntr,0)<<"=\n";
                cout<<dN_by_dEps[ElementType][GaussPntr];
            }
        }
    }
};

#endif // TRIALFUNCTION_HPP
