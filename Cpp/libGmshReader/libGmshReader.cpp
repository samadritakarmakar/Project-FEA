#ifndef LIBGMSHREADER_M_HPP
#define LIBGMSHREADER_M_HPP
#include "libGmshReader.h"



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
        int elementType[elementTypes.size()],dim2,TestStatement;
        NumOfElementTypes=elementTypes.size();
        AllocateElementData();
        std::vector<double> parametricCoord;
        //int ElementTagSize=0;
        for (int i=0;i<elementTags.size();i++)
        {
            //ElementTagSize=elementTags[i].size();
            GmshElementType[i]=elementTypes[i];
            gmsh::model::mesh::getElementProperties(GmshElementType[i], GmshElementName[i], dim2, order[i], NumOfElementNodes[i], parametricCoord);
            ElementTag[i].set_size(1,elementTags[i].size());
            GmshNodeTag[i].set_size(elementTags[i].size(),NumOfElementNodes[i]);
            cout<<"Number of Nodes per Element of Element Type "<<i+1<<" = "<<NumOfElementNodes[i]<<"\n";
            cout<<"Element Type"<<i+1<<"= "<<GmshElementName[i]<<"\n";
            cout<<"Number of Element of Tag"<<i+1<<"= "<<elementTags[i].size()<<"\n\n";
            //cout<<"ElementTags size total "<<ElementTagSize<<"\n";
        }
        cout<<"Num Of Element Types "<<NumOfElementTypes<<"\n\n";


        //ElementTag.set_size(1,elementTags.size()*elementTags[0].size());
        //GmshNodeTag.set_size(elementTags.size()*elementTags[0].size(),NumOfElementNodes);

        int m=0,n=0,c1=0;
        for (int i=0;i<elementTags.size();i++)
        {
            for(int j=0;j<elementTags[i].size();j++)
            {
               //ElementTag(0,(i+1)*j)= elementTags[i][j];
                //cout<<"c1= "<<c1;
                ElementTag[i](0,j)= elementTags[i][j];
                ++c1;
            }
            m=0;
            n=0;
            for(int k=0; k<nodeTags[i].size();k++)
            {
                //cout<<m<<","<<n<<" ";
                GmshNodeTag[i](m,n)=nodeTags[i][k];
                n++;
                TestStatement=((k+1)%NumOfElementNodes[i]==0);
                    m=m+TestStatement*1; //If NumOfElementNodes divides k+1, then increase m
                    n=!TestStatement*n;  //If NumOfElementNodes divides k+1, then increase n=0
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
        for(int j=0; j<NumOfElementTypes; ++j)
        {
            //Size of ElementNodes is same as GmshNodeTag variable
            ElementNodes[j].set_size(GmshNodeTag[j].n_rows,GmshNodeTag[j].n_cols);
            //std::cout<<"ElementType= "<<GmshElementType[j]<<"Rows= "<<GmshNodeTag[j].n_rows<<"\nColumns= "<<GmshNodeTag[j].n_cols<<"\n";
            //Arranges the unique Node tags in an assending manner
            uvec ContainsNodeTags=unique(GmshNodeTag[j]);
            ContainsNodes[j]=ContainsNodeTags;
            cout<<"Unique Nodes of Type "<<GmshElementName[j]<<" = "<<ContainsNodeTags.n_rows<<"\n";
            //uvec NodeTagPos(1);
            std::thread *MeshReaderThreads = nullptr;
            unsigned int NumofThreads=MeshReaderThreads->hardware_concurrency();
            MeshReaderThreads = new std::thread [NumofThreads];
            std::cout<<"Number of Threads to launch: "<<NumofThreads<<"\n";
            //FillElementNodes(0, ContainsNodeTags.n_rows, j, ContainsNodeTags);
            int k, stop=ContainsNodeTags.n_rows/NumofThreads;
            for (k=0; k<NumofThreads; k++)
            {
                MeshReaderThreads[k]=std::thread (&MeshReader::FillElementNodes, this, k*stop, (k+1)*stop, j, std::ref(ContainsNodeTags) );
            }
            if (ContainsNodeTags.n_rows%NumofThreads!=0)
            {
                FillElementNodes(k*stop, ContainsNodeTags.n_rows, j, ContainsNodeTags);
            }
            for (k=0; k<NumofThreads; k++)
            {
                MeshReaderThreads[k].join();
            }
            delete [] MeshReaderThreads;
            //FillElementNodes(0, ContainsNodeTags.n_rows, j, ContainsNodeTags);
            /*for (int i=0; i<(ContainsNodeTags.n_rows); i++)
            {
                //The goal is to find the same Node Tags in NodeTag and in GmshNodeTags.
                uvec GmshElemNodeTagPos=find(GmshNodeTag[j]==ContainsNodeTags(i));
                uvec NodeTagPos=find(NodeTag==ContainsNodeTags(i));
                umat temp=NodeTagPos(0)*ones<umat>(GmshElemNodeTagPos.size());
                //Here the ElementNodes are set according to index position of
                //the Node tag in the variable NodeTag. The hope is faster access
                //during solving the FEM.
                ElementNodes[j].elem(GmshElemNodeTagPos)=temp;
            }*/

        }
    }
}

void libGmshReader::MeshReader::FillElementNodes(int start, int end, int ElementType, uvec &ContainsNodeTags)
{
    for (int i=start; i<(end); i++)
    {
        //The goal is to find the same Node Tags in NodeTag and in GmshNodeTags.
        uvec GmshElemNodeTagPos=find(GmshNodeTag[ElementType]==ContainsNodeTags(i));
        uvec NodeTagPos=find(NodeTag==ContainsNodeTags(i));
        umat temp=NodeTagPos(0)*ones<umat>(GmshElemNodeTagPos.size());
        //Here the ElementNodes are set according to index position of
        //the Node tag in the variable NodeTag. The hope is faster access
        //during solving the FEM.
        ElementNodes[ElementType].elem(GmshElemNodeTagPos)=temp;
        /*for(int k=0;k<GmshElemNodeTagPos.n_rows;k++)
        {
            //Here the ElementNodes are set according to index position of
            //the Node tag in the variable NodeTag. The hope is faster access
            //during solving the FEM.
            std::cout<<"k="<<k;
            ElementNodes[ElementType](GmshElemNodeTagPos(k))=NodeTagPos(0);
        }*/
    }
}

void libGmshReader::MeshReader::setFileName(std::string FileName)
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
       throw;
    }
}

void libGmshReader::MeshReader::setDimension(int dimension)
{
    NodeData::dim=dimension;
    ElementData::dim=NodeData::dim;
}

void libGmshReader::MeshReader::FindMaxNodeNumber()
{
    maxNodeNumber=0;
    for (int i=0; i<NumOfElementTypes; i++)
    {
        int LargestNode=ContainsNodes[i](ContainsNodes[i].n_rows-1);
        if (LargestNode>maxNodeNumber)
        {
            maxNodeNumber=LargestNode;
        }
    }
    std::cout<<"Largest Node Number is = "<<maxNodeNumber<<"\n";
}

#endif
