#ifndef LIBGMSHREADER_H
#define LIBGMSHREADER_H
#include <string>
#include <armadillo>
#include <gmsh.h>
#include <string>
#include <vector>
#include <fstream>
#include <chrono>
#include <thread>
using namespace arma;
namespace libGmshReader
{

class ElementData
{
    public:
    //std::string *Type, *Degree,
    std::string *GmshElementName;
    int *NumOfElementNodes, *order, *NumOfDimension, *GmshElementType, NumOfElementTypes, maxNodeNumber;
    umat *ElementNodes, *ElementTag, *ContainsNodes, *GmshNodeTag;
    std::string fileName;
    bool fileExist;
    int dim;
    void GetElementData();
    void AllocateElementData()
    {
        //Type=new std::string [NumOfElementTypes];
        //Degree=new std::string [NumOfElementTypes];
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
        //delete []Type;
        //delete []Degree;
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
/*  ///Sets only dimension. Need to use setFileName in this case
    MeshReader(int dimension)
    {
        setDimension(dimension);
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
        FindMaxNodeNumber();
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


    ///sets file name in NodeData and ElementData
    void setFileName(std::string FileName);

    ///sets dimension in NodeData and ElementData
    void setDimension(int dimension);

    ///sets the Data of ElementNodes.
    void setElementNodes();
    ///finds max node number
    void FindMaxNodeNumber();

    int success;
//private:
    ///Fills the variable ElementNode from start to end-1.
    void FillElementNodes(int start, int end, int ElementType, uvec &ContainsNodeTags);
};
}
#endif // LIBGMSHREADER_H
