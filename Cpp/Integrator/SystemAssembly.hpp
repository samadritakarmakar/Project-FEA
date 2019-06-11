#ifndef SYSTEMASSEMBLY_HPP
#define SYSTEMASSEMBLY_HPP
#include "LocalIntegration.hpp"
#include "FEMtools.h"
template<class GenericLocalIntegrator, class GenericTrialFunction>
class SystemAssembler
{
public:
    /// Constructor. Sets the internal values of a, u &v. Also Initializes a pointer instance of Local Integrator class.
    SystemAssembler(Form<GenericTrialFunction>& a, GenericTrialFunction& u, TestFunctionGalerkin<GenericTrialFunction>& v):
        a_Internal(a), u_Internal(u), v_Internal(v)
    {
        //integrate=std::
        //        shared_ptr<LocalIntegrator<GenericTrialFunction>>
         //                                        (new LocalIntegrator<GenericTrialFunction>(a_Internal,u_Internal,v_Internal));
    }

    /// This function causes the weak form defined in the class of Integrate to be executed.
    void SetLocalIntegrator(GenericLocalIntegrator& Integrate)
    {
        integrate=&Integrate;
    }

    /// This Function is responsible for Element level Integration. At every run it Integrates one element.
    /// In the next run, it integrates the one next in the list.
    void RunLocalIntegration()
    {
        integrate->local_intergrator();
    }

    void RunLocalIntegrationVector()
    {
        integrate->local_intergrator_vector();
    }

    /// This set the matrix size of the matrix A. It is determines from the size of Function v and Function u.
    void SetMatrixSize(sp_mat& A)
    {
        A_Internal=&A;
        A_Internal->resize
                (v_Internal.vectorLvl*u_Internal.Msh->NodeTag.n_rows,u_Internal.vectorLvl*u_Internal.Msh->NodeTag.n_rows);
        //cout<<"No of Rows of NodeTags= "<<u_Internal.Msh->NodalCoordinates;
    }

    /// This set the vector size of the matrix b. It is determines from the size of Function v.
    void SetVectorSize(mat &b)
    {
        b_Internal=&b;
        b_Internal->zeros(v_Internal.vectorLvl*u_Internal.Msh->NodeTag.n_rows,1);
        std::cout<<"Size of vector b= "<<v_Internal.vectorLvl*u_Internal.Msh->NodeTag.n_rows<<", 1"<<"\n";
    }

    /// This Function sets the weak form to be inegrated at the local level. It integrated over every element
    /// and then fills up the global matrix A.
    void RunSystemAssembly(GenericLocalIntegrator& Integrate, sp_mat& A)
    {
        SetLocalIntegrator(Integrate);
        cout<<"Size of A= "<<A.n_rows<<","<<A.n_cols<<"\n";
        //Configuration for Batch addition of matrix to global matrix
        bool add_values=false;
        bool sort_locations=true;
        bool check_for_zeros=true;
        //------------------------------------------------------------
        NodePositions=std::vector<umat>(u_Internal.NoOfElementTypes);
        for (int ElmntTyp=0; ElmntTyp<u_Internal.NoOfElementTypes; ElmntTyp++)
        {
            GetNodePostions(NodePositions[ElmntTyp], u_Internal.ElmntNodes[ElmntTyp], u_Internal.vectorLvl);
            //u_Internal.Msh->ElementNodes[ElmntTyp].n_rows;
            for (int ElmntNmbr=0; ElmntNmbr<u_Internal.NoOfElements[ElmntTyp] ; ElmntNmbr++)
            {
                RunLocalIntegration();
                umat positions=NodePositions[ElmntTyp].row(a_Internal.ElementNumber);
                /*umat locations(2,positions.n_cols*positions.n_cols);
                int locationPtr=0;
                for (int row=0;row<positions.n_cols;row++)
                {
                    for (int col=0;col<positions.n_cols;col++)
                    {
                        umat locations1;
                        locations1<<positions(row)<<endr<<positions(col)<<endr;
                        locations.col(locationPtr)=locations1;
                        locationPtr++;
                    }
                }
                vec values=vectorise(mat(a_Internal.ResultingMat.t()));
                //cout<<"Size of locations ="<<locations.n_rows<<","<<locations.n_cols<<"\n";
                //cout<<"Size of values ="<<values.n_rows<<","<<values.n_cols<<"\n";
                sp_mat A_temp=sp_mat(add_values, locations, values, A.n_rows, A.n_cols, sort_locations, check_for_zeros);*/

                sp_mat A_temp=BatchFill_Atemp(A, ElmntTyp, positions);
                A=A+A_temp;

                a_Internal.NextElementNumber();
            }
            a_Internal.NextElementType();
        }
    }

    void RunSystemAssemblyVector(GenericLocalIntegrator& Integrate, mat& b)
    {
        SetLocalIntegrator(Integrate);
        NodePositions=std::vector<umat>(u_Internal.NoOfElementTypes);
        for (int ElmntTyp=0; ElmntTyp<u_Internal.NoOfElementTypes; ElmntTyp++)
        {
            GetNodePostions(NodePositions[ElmntTyp], u_Internal.ElmntNodes[ElmntTyp], u_Internal.originalVctrLvl);
            //cout<<"u_Internal.ElmntNodes[ElmntTyp] "<<u_Internal.ElmntNodes[ElmntTyp];
            //cout<<"NodePositions[ElmntTyp] "<<NodePositions[ElmntTyp]<<"\n";
            for (int ElmntNmbr=0; ElmntNmbr<u_Internal.NoOfElements[ElmntTyp]; ElmntNmbr++)
            {
                RunLocalIntegrationVector();
                umat positions=NodePositions[ElmntTyp].row(a_Internal.ElementNumber);
                umat positions2={0};
                b.submat(positions, positions2)=b.submat(positions, positions2)+a_Internal.ResultingVector;
                a_Internal.NextElementNumber();
            }
            a_Internal.NextElementType();
        }
    }


private:
    Form<GenericTrialFunction>& a_Internal;
    GenericTrialFunction& u_Internal;
    TestFunctionGalerkin<GenericTrialFunction>& v_Internal;
    LocalIntegrator<GenericTrialFunction>* integrate;
    std::vector<umat> NodePositions;
    sp_mat* A_Internal;
    mat * b_Internal;

    inline sp_mat BatchFill_Atemp(sp_mat& A, int& ElmntTyp, const umat& positions)
    {
        //Configuration for Batch addition of matrix to global matrix
        bool add_values=false;
        bool sort_locations=true;
        bool check_for_zeros=true;
        //------------------------------------------------------------
        umat locations(2,positions.n_cols*positions.n_cols);
        int locationPtr=0;
        for (int row=0;row<positions.n_cols;row++)
        {
            for (int col=0;col<positions.n_cols;col++)
            {
                umat locations1;
                locations1<<positions(row)<<endr<<positions(col)<<endr;
                locations.col(locationPtr)=locations1;
                locationPtr++;
            }
        }
        vec values=vectorise(mat(a_Internal.ResultingMat.t()));
        //cout<<mat(a_Internal.ResultingMat.t())<<"\n";
        //cout<<"Size of locations ="<<locations.n_rows<<", "<<locations.n_cols<<"\n";
        //cout<<"Size of values ="<<values.n_rows<<", "<<values.n_cols<<"\n";
        sp_mat A_temp=sp_mat(add_values, locations, values, A.n_rows, A.n_cols, sort_locations, check_for_zeros);
        return A_temp;
    }
};

#endif // SYSTEMASSEMBLY_HPP
