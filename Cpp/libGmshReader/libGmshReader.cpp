#ifndef LIBGMSHREADER_M_HPP
#define LIBGMSHREADER_M_HPP
#include "libGmshReader.h"
#include <vector>
#include <armadillo>
#include <gmsh.h>
#include <string>
#include <vector>
#include <fstream>

void libGmshReader::NodeData::GetNodeData()
{
    gmsh::initialize();
    //std::ifstream fileExist(fileName);
    if(fileExist)
    {
        gmsh::open(fileName);
        std::vector<int> nodeTags;
        std::vector<double> coord, parametricCoord;
        const int tag=-1, dim2=-1;//Retrieves all tags and dimensions
        const bool includeBoundary = true;
        gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, dim2, tag, includeBoundary);
        std::vector<int>::const_iterator j=nodeTags.begin();
        int k=1,m=0,n=0;
        NumOfNodes=nodeTags.size();
        NodeTag.set_size(NumOfNodes,1);
        NodalCoordinates.set_size(NumOfNodes,3);
        int kMod3Bool;
        for(std::vector<double>::const_iterator i=coord.begin(); i != coord.end();i++)
        {

            NodeTag(m,0)=*j;
            NodalCoordinates(m,n)=*i;
            kMod3Bool=(k%3==0);
            j=j+(kMod3Bool)*1; //Goes to next nodeTag only if 3 data points have been moved over
            m=m+(kMod3Bool)*1; //Goes to next nodeTag only if 3 data points have been moved over
            n=(++n)*!kMod3Bool; //Sets n to 0 if all 3 dimensions have been retrived
            k++;
        }
        //return 1;
    }
    //else
      //  return 0;
}

void libGmshReader::ElementData::GetElementData()
{
    if (fileExist)
    {
        std::vector<int> elementTypes;
        std::vector<std::vector<int> > elementTags, nodeTags;
        const int tag = -1;
        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, dim, tag);
        int elementType=elementTypes[0],dim2,TestStatement;
        std::vector<double> parametricCoord;
        gmsh::model::mesh::getElementProperties(elementType, GmshElementName, dim2, order, NumOfElementNodes, parametricCoord);
        ElementTag.set_size(1,elementTags.size()*elementTags[0].size());
        GmshNodeTag.set_size(elementTags.size()*elementTags[0].size(),NumOfElementNodes);
        int m=0,n=0;
        for (int i=0;i<elementTags.size();i++)
        {
            for(int j=0;j<elementTags[i].size();j++)
            {
               ElementTag(0,(i+1)*j)= elementTags[i][j];
            }
            for(int k=0; k<nodeTags[i].size();k++)
            {
                GmshNodeTag(m*(i+1),n)=nodeTags[i][k];
                n++;
                TestStatement=((k+1)%NumOfElementNodes==0);
                    m=m+TestStatement*1;
                    n=!TestStatement*n;
            }
        }
        //return 1;
    }
    /*else
    {
        return 0;
    }*/
}
void libGmshReader::MeshReader::setElementNodes()
{
    if (success)
    {
        //Size of ElementNodes is same as GmshNodeTag variable
        ElementNodes.set_size(GmshNodeTag.n_rows,GmshNodeTag.n_cols);
        //Arranges the unique Node tags in an assending manner
        uvec ContainsNodeTags=unique(GmshNodeTag);
        uvec NodeTagPos(1);
        for (int i=0; i<(ContainsNodeTags.n_rows); i++)
        {
            //The goal is to find the same Node Tags in NodeTag and in GmshNodeTags.
            uvec GmshElemNodeTagPos=find(GmshNodeTag==ContainsNodeTags(i));
            NodeTagPos=find(NodeTag==ContainsNodeTags(i));
            umat temp=NodeTagPos(0)*ones<umat>(GmshElemNodeTagPos.size());
            ElementNodes.elem(GmshElemNodeTagPos)=temp;
            /*for(int j=0;j<GmshElemNodeTagPos.n_rows;j++)
            {
                //Here the ElementNodes are set according to index position of
                //the Node tag in the variable NodeTag. The hope is faster access
                //during solving the FEM.
                ElementNodes(GmshElemNodeTagPos(j))=NodeTagPos(0);
            }*/
        }
    }

}
#endif
