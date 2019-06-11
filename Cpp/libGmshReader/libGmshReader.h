#ifndef LIBGMSHREADER_H
#define LIBGMSHREADER_H
#include <string>
#include <armadillo>
#include "gmsh.h"
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
    std::vector<std::string> GmshElementName, GmshElementNameOnly;
    int  NumOfElementTypes, maxNodeNumber, NumOfPhysclGrps;
    std::vector<int> NumOfElementNodes, order, NumOfDimension, GmshElementType;
    std::vector <umat> ElementNodes, ElementTag, ContainsNodes, GmshNodeTag;
    std::vector <std::vector <umat>> ElmntPhysclGrpNodes;
    std::vector<std::string> PhysicalGroupName;
    std::vector<size_t> PhysicalGroupNodeTags;
    std::string fileName;
    bool fileExist;
    int dim;
protected:
    void GetElementData();
    /// Allocate the element data as per the number of element types.
    void AllocateElementData()
    {
        GmshElementName =std::vector <std::string> (NumOfElementTypes);
        GmshElementNameOnly =std::vector <std::string> (NumOfElementTypes);
        NumOfElementNodes =std::vector <int> (NumOfElementTypes);
        order=std::vector <int> (NumOfElementTypes);
        NumOfDimension=std::vector <int> (NumOfElementTypes);
        GmshElementType=std::vector <int> (NumOfElementTypes);
        ElementNodes=std::vector <arma::umat> (NumOfElementTypes);
        ElementTag=std::vector <arma::umat> (NumOfElementTypes);
        ContainsNodes=std::vector <arma::umat> (NumOfElementTypes);
        GmshNodeTag=std::vector <arma::umat> (NumOfElementTypes);
        ElmntPhysclGrpNodes =std::vector <std::vector <umat>> (NumOfElementTypes);
    }
    /// Delete Element data
    void DeleteElementData()
    {
        //gmsh::finalize();
    }
    /// Extract just the element name and remove the number of nodes from it.
    void GetGmshElementNameOnly();
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
protected:
    ///Extract node Data
    void GetNodeData();

};
class MeshReader: public ElementData, public NodeData
{
public:

   ///Need to use setFileName in this case
    MeshReader()
    {
    }
 /*    ///Sets only fileName. Need to use setDimension in this case.
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
        /// Extract just the element name and remove the number of nodes from it.
        GetGmshElementNameOnly();
        ///Sets the variable ElementNodes Mesh file
        setElementNodes();
        FindMaxNodeNumber();
        end=std::chrono::system_clock::now();
        std::cout<<"Done Reading the Mesh!\n";
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout<<"Time taken to Read Mesh= "<<elapsed_seconds.count()<<" seconds\n";
        GetPhysicalGroupData();
    }
    /// This Copies the FileName, Node Data from the given File and skips reading Node Data
    MeshReader(libGmshReader::MeshReader GivenMsh, int dimension)
    {
        std::cout<<"Reading the Mesh...\n";
        ///Sets fileName
        setFileName(GivenMsh.NodeData::fileName);
        ///sets dimension
        setDimension(dimension);
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        gmsh::initialize();
        NumOfNodes=GivenMsh.NumOfNodes;
        NodeTag.set_size(GivenMsh.NodeTag.n_rows, GivenMsh.NodeTag.n_cols);
        NodeTag=GivenMsh.NodeTag;
        NodalCoordinates.set_size(GivenMsh.NodalCoordinates.n_rows,GivenMsh.NodalCoordinates.n_cols);
        NodalCoordinates=GivenMsh.NodalCoordinates;
        gmsh::open(ElementData::fileName);
        ///Extracts Element data from Mesh file
        GetElementData();
        /// Extract just the element name and remove the number of nodes from it.
        GetGmshElementNameOnly();
        ///Sets the variable ElementNodes Mesh file
        setElementNodes();
        FindMaxNodeNumber();
        end=std::chrono::system_clock::now();
        std::cout<<"Done Reading the Mesh!\n";
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout<<"Time taken to Read Mesh= "<<elapsed_seconds.count()<<" seconds\n";
        GetPhysicalGroupData();
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

    void GetPhysicalGroupData();

    int success;


private:
    ///Fills the variable ElementNode from start to end-1.
    void FillElementNodes(int start, int end, int ElementType, uvec &ContainsNodeTags);
    gmsh::vectorpair dimTags;

};
}
#endif // LIBGMSHREADER_H
