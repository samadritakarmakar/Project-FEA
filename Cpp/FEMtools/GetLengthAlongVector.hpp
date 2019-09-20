#ifndef GETLENGTHALONGVECTOR_HPP
#define GETLENGTHALONGVECTOR_HPP

#include "../Math/vecnorm.hpp"
//#include "../TestAndTrialFunction/TrialFunction.hpp"
using namespace arma;
template<class GenericTrialFunction>
mat GetLengthAlongVector(GenericTrialFunction &u, const int Dimension, const vec &vector)
{
    int ElmntNmbr=0, ElmntType=0;
    vec vel=vector.rows(0,Dimension-1);
    mat Vector, max_cos_theta(1,1);
    mat cos_theta;
    for (ElmntType=0;ElmntType<u.NoOfElementTypes; ElmntType++)
    {
        Vector.resize(ElmntType*u.NoOfElements[ElmntType]+u.NoOfElements[ElmntType],Dimension);
        for (ElmntNmbr=0; ElmntNmbr<u.NoOfElements[ElmntType]; ElmntNmbr++)
        {
            umat NodesAtElmnt = u.GetNodesAt(0, ElmntNmbr);
            mat Coordinates1 = u.GetCoordinatesAt(NodesAtElmnt);
            mat Coordinates =Coordinates1.cols(0, Dimension-1);
            mat VectorTransposed(Coordinates.n_rows, Coordinates.n_cols);
            VectorTransposed.rows(1, Coordinates.n_rows-1) = Coordinates.rows(1, Coordinates.n_rows-1)-Coordinates.rows(0, Coordinates.n_rows-2);
            VectorTransposed.row(0) = Coordinates.row(0)-Coordinates.row(Coordinates.n_rows-1);
            //cout<<"VectorTransposed= "<<VectorTransposed<<"\n";
            vec Vector_dot_Vel=VectorTransposed*vel;
            cos_theta.set_size(Vector_dot_Vel.n_rows, Vector_dot_Vel.n_cols);
            cos_theta=Vector_dot_Vel/(vecnorm(VectorTransposed.t())*norm(vel));
            //cout<<"cos_theta = "<<cos_theta<<"\n";
            max_cos_theta=max(cos_theta);
            umat index_maxm=index_max(abs(cos_theta));
            Vector.row(ElmntType*ElmntNmbr+ElmntNmbr)=(VectorTransposed.rows(index_maxm));
        }
    }
    return Vector;
}


#endif
