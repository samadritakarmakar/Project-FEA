#ifndef LIBGMSHREADER_H
#define LIBGMSHREADER_H
#include <string>
#include <armadillo>
#include <gmsh.h>
#include <string>
#include <vector>
#include <fstream>
#include <chrono>
using namespace arma;
namespace libGmshReader
{

class ElementData
{
    public:
    std::string *Type, *Degree, *GmshElementName;
    int *NumOfElementNodes, *order, *NumOfDimension, *GmshElementType, NumOfElementTypes;
    umat *ElementNodes, *ElementTag, *ContainsNodes, *GmshNodeTag;
    std::string fileName;
    bool fileExist;
    int dim;
    void GetElementData();
    void AllocateElementData()
    {
        Type=new std::string [NumOfElementTypes];
        Degree=new std::string [NumOfElementTypes];
        GmshElementName=new std::string [NumOfElementTypes];
        NumOfElementNodes =new int [NumOfElementTypes];
        order=new int [NumOfElementTypes];
        NumOfDimension=new int [NumOfElementTypes];
        GmshElementType=new int [NumOfElementTypes];
        ElementNodes=new umat [NumOfElementTypes];
        ElementTag=new umat [NumOfElementTypes];
        ContainsNodes=new umat [NumOfElementTypes];
        GmshNodeTag=new umat[NumOfElementTypes];
    }
    void DeleteElementData()
    {
        delete []Type;
        delete []Degree;
        delete []GmshElementName;
        delete []NumOfElementNodes;
        delete []order;
        delete []NumOfDimension;
        delete []ElementNodes;
        delete []ElementTag;
        delete []ContainsNodes;
        delete []GmshNodeTag;
    }

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
    void GetNodeData();

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
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        ///Extracts Node data from Mesh file
        GetNodeData();
        ///Extracts Element data from Mesh file
        GetElementData();
        ///Sets the variable ElementNodes Mesh file
        setElementNodes();
        end=std::chrono::system_clock::now();
        std::cout<<"Done Reading the Mesh!\n";
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout<<"Time taken to Read Mesh= "<<elapsed_seconds.count()<<" seconds\n";
    }
    ///Deallocates all allocated Element Data
    ~MeshReader()
    {
        DeleteElementData();
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
