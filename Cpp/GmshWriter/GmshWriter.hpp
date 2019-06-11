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
    GmshWriter(TrialFunction& u, std::string OutputFileName):
        u_Internal(u), outputFileName(OutputFileName)
    {
        modelName=GetOnlyFileName();
    }

    void GoToNextViewTag()
    {
        currentViewTag++;
    }
    void WriteToGmsh(mat &x, double time=0.0)
    {
        gmsh::initialize();
        gmsh::open(u_Internal.Msh->ElementData::fileName);
        gmsh::view::add(modelName, currentViewTag);
        std::string dataType="NodeData";
        std::vector<std::size_t> Nodetags;
        std::vector<std::vector<double>> Nodedata;
        getNodeTagAndNodeData(x, Nodetags, Nodedata);
        gmsh::view::addModelData(currentViewTag, 0, modelName, dataType, Nodetags, Nodedata, time);
        gmsh::view::write(currentViewTag, outputFileName);
    }
private:
    TrialFunction& u_Internal;
    std::string& outputFileName;
    std::string modelName;
    int currentViewTag=1;

    std::string GetOnlyFileName()
    {
        std::string FileName=u_Internal.Msh->NodeData::fileName;
        int EndOfSlash=FileName.find_last_of("/");
        int EndOfDot=FileName.find_last_of(".");
        return FileName.substr(EndOfSlash+1, EndOfDot-EndOfSlash-1);
    }

    void getNodeTagAndNodeData(mat &x, std::vector<std::size_t>& Nodetags, std::vector<std::vector<double>>& Nodedata)
    {
        umat& NodeTagPtr= u_Internal.Msh->NodeTag;
        Nodetags= std::vector<std::size_t> (NodeTagPtr.n_rows);
        Nodedata= std::vector<std::vector<double>>  (NodeTagPtr.n_rows);
        int& u_originalVctrLvl=u_Internal.originalVctrLvl;
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
