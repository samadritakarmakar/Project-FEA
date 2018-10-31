#ifndef LIBMESHREADER_H
#define LIBMESHREADER_H
#include <armadillo>
#include <string>
using namespace arma;
using namespace std;
class ElementType
{
public:
    string Type, Degree;
    int NumOfElementNodes, NumofDimensions, ContainsNodes;
    mat ElementNodes;
};
class MeshReader : public ElementType
{
  public:

    void GetElemProperty(int GmshElementNum, string& Type, string& Degree, int& NumOfElementNodes, int& NumofDimensions, int& Supported);
    mat SortAndKeep(mat ElementNodes);

};
//mat SortAndKeep(mat ElementNodes);
#endif // LIBMESHREADER_H
