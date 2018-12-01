#ifndef LIBGMSHREADER_H
#define LIBGMSHREADER_H
#include <string>
#include <armadillo>
#include <gmsh.h>
#include <string>
#include <vector>
#include <fstream>
using namespace arma;
namespace libGmshReader
{

class ElementData
{
    public:
    std::string Type, Degree, GmshElementName;
    int NumOfElementNodes, order, NumOfDimension;
    umat ElementNodes, ElementTag, ContainsNodes, GmshNodeTag;
    std::string fileName;
    bool fileExist;
    int dim;
    int GetElementData();
};
class NodeData
{
public:
    int NumOfNodes;
    umat NodeTag;
    mat NodalCoordinates;
    std::string fileName;
    bool fileExist;
    int dim;
    int GetNodeData();

};
class MeshReader: public ElementData, public NodeData
{
public:

/*    ///Need to use setFileName in this case
    MeshReader()
    {
    }
    ///Sets only fileName. Need to use setDimension in this case.
    MeshReader(std::string FileName)
    {
        setFileName(FileName);
    }*/
    ///Sets fileName and Dimension.
    MeshReader(std::string FileName, int dimension)
    {
        std::cout<<"Reading the Mesh...\n";
        ///Sets fileName
        setFileName(FileName);
        ///sets dimension
        setDimension(dimension);
        ///Extracts Node data from Mesh file
        GetNodeData();
        ///Extracts Element data from Mesh file
        GetElementData();
        ///Sets the variable ElementNodes Mesh file
        setElementNodes();
        std::cout<<"Done Reading the Mesh!\n";
    }
  /*  ///Sets only dimension. Need to use setFileName in this case
    MeshReader(int dimension)
    {
        setDimension(dimension);
    }*/

    ///sets file name in NodeData and ElementData
    void setFileName(std::string FileName)
    {
        ///Check if file Exists.
        std::ifstream fileExists(FileName);
        if(fileExists)
        {
           NodeData::fileName=FileName;
           NodeData::fileExist=true;
           ElementData::fileName=NodeData::fileName;
           ElementData::fileExist=NodeData::fileExist;
           success=1;
        }
        else
        {
           std::cout<<"File Does not exist!\n";
           NodeData::fileExist=false;
           ElementData::fileExist=NodeData::fileExist;
           success=0;
        }
    }
    ///sets dimension in NodeData and ElementData
    void setDimension(int dimension)
    {
        NodeData::dim=dimension;
        ElementData::dim=NodeData::dim;
    }
    ///sets the Data of ElementNodes.
    void setElementNodes();

    int success;
};
}
#endif // LIBGMSHREADER_H
