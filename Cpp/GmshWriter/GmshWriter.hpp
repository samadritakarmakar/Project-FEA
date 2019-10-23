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
    std::string dataType="NodeData";
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
        //std::string dataType="NodeData";
        std::vector<std::size_t> tags;
        std::vector<std::vector<double>> data;
        int numberOfTags;
        if(dataType.compare("NodeData")==0)
        {
            numberOfTags=GetNodeTags();
        }
        else if(dataType.compare("ElementData")==0)
        {
            numberOfTags=GetElementTags();
        }
        if(dataType.compare("NodeData")==0)
        {
            getNodeTagAndNodeData(x, tags, data);
        }
        else if(dataType.compare("ElementData")==0)
        {
            getElementTagAndElementData(x, tags, data);
        }
        gmsh::view::addModelData(currentViewTag, step, modelName, dataType, tags, data, time, u_Internal.originalVctrLvl);
        if (step==0)
        {
            gmsh::view::write(currentViewTag, outputFileName);
        }
        else
        {
            gmsh::view::write(currentViewTag, outputFileName, true);
        }

        gmsh::finalize();
    }

    void WriteToGmshSymmetric_3x3(mat &x, int step=0, double time=0.0)
    {
        gmsh::initialize();
        gmsh::open(u_Internal.Msh->ElementData::fileName);
        //gmsh::view::add(modelName, currentViewTag);
        gmsh::view::add(viewName, currentViewTag);
        //std::string dataType="NodeData";
        std::vector<std::size_t> tags;
        std::vector<std::vector<double>> data;
        mat I=eye(6,6);
        umat Positions={{0,3,5,3,1,4,5,4,2}};
        //This is a permutaion matrix. This is used to make a symmetric (voigt notation) with 6 components back
        //to a 9 component representation.
        sp_mat PermMat=sp_mat(I.cols(Positions));
        int numberOfTags;
        if(dataType.compare("NodeData")==0)
        {
            numberOfTags=GetNodeTags();
        }
        else if(dataType.compare("ElementData")==0)
        {
            numberOfTags=GetElementTags();
        }
        mat x_new(numberOfTags*9,1);
        for (int node = 0; node <numberOfTags ; ++node)
        {
           x_new(span(node*9,node*9+8),0)= PermMat.t()*x(span(node*6,node*6+5),0);
        }
        if(dataType.compare("NodeData")==0)
        {
            getNodeTagAndNodeData(x_new, tags, data, 9);
        }
        else if(dataType.compare("ElementData")==0)
        {
            getElementTagAndElementData(x_new, tags, data, 9);
        }
        gmsh::view::addModelData(currentViewTag, step, modelName, dataType, tags, data, time, 9);
        gmsh::view::write(currentViewTag, outputFileName);
        gmsh::finalize();
    }

    void SetDataType_to_NodeData()
    {
        dataType="NodeData";
    }

    void SetDataType_to_ElementData()
    {
        dataType="ElementData";
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

    /// Sets NodeData and NodeTags and returns the total number of Node tags (Total Number of Nodes).
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
        //cout<<"NodeTagPtr.n_rows ="<<NodeTagPtr.n_rows<<"\n";
        for (int tagcount=0; tagcount<NodeTagPtr.n_rows; tagcount++)
        {
            Nodetags[tagcount]=NodeTagPtr(tagcount,0);
            Nodedata[tagcount]=std::vector<double>(u_originalVctrLvl);
            for (int vctrLvlCntr=0; vctrLvlCntr<u_originalVctrLvl; vctrLvlCntr++)
            {
                //cout<<" at [tagcount]= "<<tagcount<<" [vctrLvlCntr]= "<<vctrLvlCntr<<"\n";
                Nodedata[tagcount][vctrLvlCntr]=x(tagcount*u_originalVctrLvl+vctrLvlCntr);
                //cout<<Nodedata[tagcount][vctrLvlCntr]<<" ";
            }

            //cout<<"\n";
        }
    }

    /// Sets ElementData and ElementTags and returns the total number of Element tags.
    void getElementTagAndElementData(mat &x, std::vector<std::size_t>& ElementTags, std::vector<std::vector<double>>& ElementData,
                               int New_vctrLvl=0)
    {
        std::vector<umat> ElementTagPtr(u_Internal.Msh->ElementTag.size());
        int ElementTagSize=0;
        for(int ElmntType=0; ElmntType<u_Internal.NoOfElementTypes; ElmntType++)
        {
            ElementTagPtr[ElmntType]=u_Internal.Msh->ElementTag[ElmntType];
            ElementTagSize=ElementTagSize+ElementTagPtr[ElmntType].n_cols;
        }
        ElementTags= std::vector<std::size_t> (ElementTagSize);
        ElementData= std::vector<std::vector<double>>  (ElementTagSize);
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
        int ElmntType=0;
        for (int tagcount=0; tagcount<ElementTagSize; tagcount++)
        {
            ElementTags[tagcount]=ElementTagPtr[ElmntType](0,tagcount);
            if(ElementTagPtr[ElmntType].n_cols-1==tagcount)
            {
                //Increases the ElmntType count by 1 when it reaches the end of the current ElementType.
                ElmntType++;
            }
            ElementData[tagcount]=std::vector<double>(u_originalVctrLvl);
            for (int vctrLvlCntr=0; vctrLvlCntr<u_originalVctrLvl; vctrLvlCntr++)
            {
                ElementData[tagcount][vctrLvlCntr]=x(tagcount*u_originalVctrLvl+vctrLvlCntr);
                //cout<<Nodetags[tagcount]<<"    "<<Nodedata[tagcount][vctrLvlCntr]<<"\n";
                //cout<<Nodedata[tagcount][vctrLvlCntr]<<" ";
            }

            //cout<<"\n";
        }
    }

    int GetNodeTags()
    {
        umat& NodeTagPtr= u_Internal.Msh->NodeTag;
        return NodeTagPtr.n_rows;
    }

    int GetElementTags()
    {
        std::vector<umat> ElementTagPtr(u_Internal.Msh->ElementTag.size());
        int ElementTagSize=0;
        for(int ElmntType=0; ElmntType<u_Internal.NoOfElementTypes; ElmntType++)
        {
            ElementTagPtr[ElmntType]=u_Internal.Msh->ElementTag[ElmntType];
            ElementTagSize=ElementTagSize+ElementTagPtr[ElmntType].n_cols;
        }
        return ElementTagSize;
    }
};

#endif // GMSHWRITER_HPP
