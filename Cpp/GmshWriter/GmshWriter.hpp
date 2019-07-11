#ifndef GMSHWRITER_HPP
#define GMSHWRITER_HPP
#include <string>
#include <armadillo>
#include <gmsh.h>
#include "TrialFunction.hpp"
using namespace arma;
class GmshWriter
{
public:
    std::string viewName="Results";
    GmshWriter(TrialFunction& u, std::string OutputFileName):
        u_Internal(u), outputFileName(OutputFileName)
    {
        modelName=GetOnlyFileName();
    }

    /// Used to switch to next View tag;
    void GoToNextViewTag()
    {
        currentViewTag++;
    }
    void WriteToGmsh(mat &x, int step=0, double time=0.0)
    {
        gmsh::initialize();
        gmsh::open(u_Internal.Msh->ElementData::fileName);
        //gmsh::view::add(modelName, currentViewTag);
        gmsh::view::add(viewName, currentViewTag);
        std::string dataType="NodeData";
        std::vector<std::size_t> Nodetags;
        std::vector<std::vector<double>> Nodedata;
        getNodeTagAndNodeData(x, Nodetags, Nodedata);
        gmsh::view::addModelData(currentViewTag, step, modelName, dataType, Nodetags, Nodedata, time, u_Internal.originalVctrLvl);
        gmsh::view::write(currentViewTag, outputFileName);
        gmsh::finalize();
    }

    void WriteToGmshSymmetric_3x3(mat &x, int step=0, double time=0.0)
    {
        gmsh::initialize();
        gmsh::open(u_Internal.Msh->ElementData::fileName);
        //gmsh::view::add(modelName, currentViewTag);
        gmsh::view::add(viewName, currentViewTag);
        std::string dataType="NodeData";
        std::vector<std::size_t> Nodetags;
        std::vector<std::vector<double>> Nodedata;
        mat I=eye(6,6);
        umat Positions={{0,3,5,3,1,4,5,4,2}};
        //This is a permutaion matrix. This is used to make a symmetric (voigt notation) with 6 components back
        //to a 9 component representation.
        sp_mat PermMat=sp_mat(I.cols(Positions));
        int numOfNodes=u_Internal.Msh->NodeTag.n_rows;
        mat x_new(numOfNodes*9,1);
        for (int node = 0; node <numOfNodes ; ++node)
        {
           x_new(span(node*9,node*9+8),0)= PermMat.t()*x(span(node*6,node*6+5),0);
        }
        getNodeTagAndNodeData(x_new, Nodetags, Nodedata, 9);
        gmsh::view::addModelData(currentViewTag, step, modelName, dataType, Nodetags, Nodedata, time, 9);
        gmsh::view::write(currentViewTag, outputFileName);
        gmsh::finalize();
    }

private:
    TrialFunction& u_Internal;
    std::string outputFileName;
    std::string modelName;
    int currentViewTag=1;

    std::string GetOnlyFileName()
    {
        std::string FileName=u_Internal.Msh->NodeData::fileName;
        int EndOfSlash=FileName.find_last_of("/");
        int EndOfDot=FileName.find_last_of(".");
        return FileName.substr(EndOfSlash+1, EndOfDot-EndOfSlash-1);
    }

    void getNodeTagAndNodeData(mat &x, std::vector<std::size_t>& Nodetags, std::vector<std::vector<double>>& Nodedata,
                               int New_vctrLvl=0)
    {
        umat& NodeTagPtr= u_Internal.Msh->NodeTag;
        Nodetags= std::vector<std::size_t> (NodeTagPtr.n_rows);
        Nodedata= std::vector<std::vector<double>>  (NodeTagPtr.n_rows);
        int u_originalVctrLvl;
        if (New_vctrLvl==0)
        {
            u_originalVctrLvl=u_Internal.originalVctrLvl;
        }
        else
        {
            u_originalVctrLvl=New_vctrLvl;
        }

        cout<<"Vector Level for write ="<<u_originalVctrLvl<<"\n";
        //cout<<"Unique Node Tags are"<<unique(NodeTagPtr);
        for (int tagcount=0; tagcount<NodeTagPtr.n_rows; tagcount++)
        {
            Nodetags[tagcount]=NodeTagPtr(tagcount,0);
            Nodedata[tagcount]=std::vector<double>(u_originalVctrLvl);
            for (int vctrLvlCntr=0; vctrLvlCntr<u_originalVctrLvl; vctrLvlCntr++)
            {
                Nodedata[tagcount][vctrLvlCntr]=x(tagcount*u_originalVctrLvl+vctrLvlCntr);
                //cout<<Nodetags[tagcount]<<"    "<<Nodedata[tagcount][vctrLvlCntr]<<"\n";
                //cout<<Nodedata[tagcount][vctrLvlCntr]<<" ";
            }

            //cout<<"\n";
        }
    }
};

#endif // GMSHWRITER_HPP
