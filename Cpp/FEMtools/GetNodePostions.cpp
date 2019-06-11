
#include "FEMtools.h"

#include <armadillo>
using namespace arma ;
/// Get Total number of degree of freedom (including vector level, 1d, 2d, 3d)
/// and accordingly set the position of the quantity in the matrix
void GetNodePostions(umat &NodePositons, const umat &ElementNodes, int vectorLevel)
{
  int column, i, j, row ;
  row = ElementNodes.n_rows ;
  column = ElementNodes.n_cols ;
  //NodePositons = arma::zeros<umat>(row, vectorLevel*column) ;
  NodePositons.set_size(row, vectorLevel*column);
  for (i=0; i<column; i++)
  {
    for (j=0; j<=vectorLevel-1; j++)
    {
       uvec oneVector= arma::ones<uvec>(row);
      //std::cout<<"Col Position= "<<vectorLevel*(i+1)-j-1<<"\n";
      //Too long formula. Any shorter formula is a welcome!
      /*Fills  the NodePostions according to vector level
        Eg: for Node Postion 0, and vectorLevel 3, values are filled at 0, 1, 2,
        for Node Postion 1, and vectorLevel 3, values are filled at 3, 4, 5.

        If at the certain ElementNodes column the node number is 0, for vectorLevel 3, then the values 0, 1, 2 are filled
        If at the certain ElementNodes column the node number is 1, for vectorLevel 3, then the values 2, 3, 4 are filled
        The general formula used is vectorLevel*(i+1)-j-1*/
      NodePositons.col(vectorLevel*(i+1)-j-1) = vectorLevel*(ElementNodes.col(i)+1*oneVector)-j*oneVector
              -1*oneVector;
    }
  }
}

void GetNodePostions(umat &NodePositons, const umat &ElementNodes, int vectorLevel, int originalVctrLvl)
{
  int column, i, j, row ;
  row = ElementNodes.n_rows ;
  column = ElementNodes.n_cols ;
  //NodePositons = arma::zeros<umat>(row, vectorLevel*column) ;
  NodePositons.set_size(row, vectorLevel*column);
  for (i=0; i<column; i++)
  {
    for (j=0; j<=vectorLevel-1; j++)
    {
       uvec oneVector= arma::ones<uvec>(row);
      //std::cout<<"Col Position= "<<vectorLevel*(i+1)-j-1<<"\n";
      //Too long formula. Any shorter formula is a welcome!
      /*Fills  the NodePostions according to vector level
        Eg: for Node Postion 0, and vectorLevel 3, values are filled at 0, 1, 2,
        for Node Postion 1, and vectorLevel 3, values are filled at 3, 4, 5.

        If at the certain ElementNodes column the node number is 0, for vectorLevel 3, then the values 0, 1, 2 are filled
        If at the certain ElementNodes column the node number is 1, for vectorLevel 3, then the values 2, 3, 4 are filled
        The general formula used is vectorLevel*(i+1)-j-1*/
      NodePositons.col(vectorLevel*(i+1)-j-1) = originalVctrLvl*(ElementNodes.col(i)+1*oneVector)-j*oneVector
              -1*oneVector;
    }
  }
}
