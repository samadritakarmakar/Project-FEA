#include "ProjectFEA.hpp"

int main(int argc, char *argv[])
{
    if(argc==1 || argc>3)
    {
        std::cout<<"Usage: ./Trial <.msh Filename> <Dimension>\n";
        return 0;
    }
    std::string FileName(argv[1]);
    int Dimension=*argv[2]-'0';
    int vectorLevel=1;
    libGmshReader::MeshReader Mesh(FileName, Dimension);
    libGmshReader::MeshReader Mesh_order_1(Mesh, Dimension, 1);
    TrialFunction u(Mesh_order_1, vectorLevel);
    vec vel1;
    vel1<<4<<endr<<3<<endr<<1;
    /*vec vel=vel1.rows(0,Dimension-1);
    for (int ElmntNmbr=0; ElmntNmbr<u.NoOfElements[0]; ElmntNmbr++)
    {
        umat NodesAtElmnt = u.GetNodesAt(0, ElmntNmbr);
        mat Coordinates1 = u.GetCoordinatesAt(NodesAtElmnt);
        mat Coordinates =Coordinates1.cols(0, Dimension-1);
        mat VectorTransposed(Coordinates.n_rows, Coordinates.n_cols);
        VectorTransposed.rows(1, Coordinates.n_rows-1) = Coordinates.rows(1, Coordinates.n_rows-1)-Coordinates.rows(0, Coordinates.n_rows-2);
        VectorTransposed.row(0) = Coordinates.row(0)-Coordinates.row(Coordinates.n_rows-1);
        //cout<<"VectorTransposed= "<<VectorTransposed<<"\n";
        vec Vector_dot_Vel=VectorTransposed*vel;
        mat cos_theta=Vector_dot_Vel/(vecnorm(VectorTransposed.t())*norm(vel));
        //cout<<"cos_theta = "<<cos_theta<<"\n";
        umat index_maxm=index_max(abs(cos_theta));
        cout<<"VectorTransposed "<<VectorTransposed.rows(index_maxm)<<"\n";
    }*/
    vec Size=GetLengthAlongVector(u,Dimension, vel1);
    cout<<Size;
}

