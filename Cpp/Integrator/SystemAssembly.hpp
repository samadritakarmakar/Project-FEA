#ifndef SYSTEMASSEMBLY_HPP
#define SYSTEMASSEMBLY_HPP
#include "LocalIntegration.hpp"
#include "FEMtools.h"
#include <thread>
template<class GenericLocalIntegrator, class GenericTrialFunction>
class SystemAssembler
{
public:
    /// Constructor. Sets the internal values of a, u &v. Also Initializes a pointer instance of Local Integrator class.
    SystemAssembler(Form<GenericTrialFunction>& a, GenericTrialFunction& u, TestFunctionGalerkin<GenericTrialFunction>& v):
        u_Internal(u), v_Internal(v)
    {
        a_Internal.SetNumOfThreads(1);
        a_Internal[0]=a;
        //integrate=std::
        //        shared_ptr<LocalIntegrator<GenericTrialFunction>>
         //                                        (new LocalIntegrator<GenericTrialFunction>(a_Internal,u_Internal,v_Internal));
    }

    /// Constructor. Sets the internal values of a, u &v. Also Initializes a pointer instance of Local Integrator class.
    SystemAssembler(FormMultiThread<GenericTrialFunction>& a, GenericTrialFunction& u,
                    TestFunctionGalerkin<GenericTrialFunction>& v):
        a_Internal(a), u_Internal(u), v_Internal(v)
    {
        numOfThreads=a_Internal.GetNumOfThreads();
        cout<<"Number of Threads to be Launched = "<<numOfThreads<<"\n";
    }


    /// This set the matrix size of the matrix A. It is determines from the size of Function v and Function u.
    void SetMatrixSize(sp_mat& A)
    {
        A_Internal=&A;
        //A_Internal->resize
        A_Internal->set_size(v_Internal.vectorLvl*v_Internal.numOfNodes,u_Internal.vectorLvl*u_Internal.Msh->NodeTag.n_rows);
                //(v_Internal.vectorLvl*u_Internal.Msh->NodeTag.n_rows,u_Internal.vectorLvl*u_Internal.Msh->NodeTag.n_rows);
        cout<<"Size of A= "<<A.n_rows<<","<<A.n_cols<<"\n";
    }

    /// This set the vector size of the matrix b. It is determines from the size of Function v.
    void SetVectorSize(mat &b)
    {
        b_Internal=&b;
        b_Internal->zeros(v_Internal.vectorLvl*v_Internal.numOfNodes, 1);
        std::cout<<"Size of vector b= "<<b.n_rows<<","<<b.n_cols<<"\n";
        //b_Internal->zeros(v_Internal.vectorLvl*u_Internal.Msh->NodeTag.n_rows,1);
        //std::cout<<"Size of vector b= "<<v_Internal.vectorLvl*u_Internal.Msh->NodeTag.n_rows<<", 1"<<"\n";
    }

    void SetScalarSize(vec &I)
    {
        int TotalNumberOfElements=0;
        for (int ElmntTyp =0; ElmntTyp<u_Internal.NoOfElementTypes;ElmntTyp++)
        {
           TotalNumberOfElements=TotalNumberOfElements+ u_Internal.NoOfElements[ElmntTyp];
        }
        cout<<"Total Number of Elements= "<<TotalNumberOfElements<<"\n";
        I.set_size(TotalNumberOfElements,1);
    }

    /// This Function sets the weak form to be inegrated at the local level. It integrated over every element
    /// and then fills up the global matrix A.
    void RunSystemAssembly(GenericLocalIntegrator& Integrate, sp_mat& A)
    {
        SetLocalIntegrator(Integrate);
        cout<<"Integrating over "<<u_Internal.PhysicalGrpName<<"\n";
        NodePositions_u=std::vector<umat>(u_Internal.NoOfElementTypes);
        NodePositions_v=std::vector<umat>(u_Internal.NoOfElementTypes);
        for (int ElmntTyp=0; ElmntTyp<u_Internal.NoOfElementTypes; ElmntTyp++)
        {
                GetNodePostions(NodePositions_u[ElmntTyp], u_Internal.ElmntNodes[ElmntTyp], u_Internal.originalVctrLvl);
                GetNodePostions(NodePositions_v[ElmntTyp], v_Internal.ElmntNodes[ElmntTyp], v_Internal.originalVctrLvl);
                //cout<<"v_Internal.originalVctrLvl ="<<v_Internal.originalVctrLvl<<"\n";
                //cout<<"u_Internal.ElmntNodes[ElmntTyp] =\n"<<u_Internal.ElmntNodes[ElmntTyp]<<"\n";
                //cout<<"NodePositions_v[ElmntTyp] =\n"<<NodePositions_v[ElmntTyp]<<"\n";
            //u_Internal.Msh->ElementNodes[ElmntTyp].n_rows;
            std::vector<sp_mat> Atemp(numOfThreads);
            int ElmntDivision=u_Internal.NoOfElements[ElmntTyp]/numOfThreads;
            if(integrate->IsUsingFormMultiThread())
            {
               Atemp=std::vector<sp_mat>(numOfThreads);
               std::thread Fill_A_Thread[numOfThreads];
               for (int thread=0; thread<numOfThreads; thread++)
               {
                   if(thread<numOfThreads-1)
                   {
                       Fill_A_Thread[thread]=std::thread(&SystemAssembler::IntegrateOverElements4Matrix, this, std::ref(A),
                                                         std::ref(Atemp[thread]), thread, ElmntTyp,
                                                         ElmntDivision*thread, ElmntDivision*(thread+1));
                   }
                   else
                   {
                       Fill_A_Thread[thread]=std::thread(&SystemAssembler::IntegrateOverElements4Matrix, this, std::ref(A),
                                                         std::ref(Atemp[thread]), thread, ElmntTyp,
                                                         ElmntDivision*thread, u_Internal.NoOfElements[ElmntTyp]);
                   }
               }
               for (int thread=0; thread<numOfThreads; thread++)
               {
                   Fill_A_Thread[thread].join();
                   A=A+Atemp[thread];
               }
            }
            else
            {
                Atemp=std::vector<sp_mat>(1);
                IntegrateElement4Matrix(A, Atemp[0], ElmntTyp);
                A=A+Atemp[0];
            }

        }
    }



    /// This is used for integrating to produce a vector. Generally used for source terms (body forces) and
    /// neumann conditions (traction forces).
    void RunSystemAssemblyVector(GenericLocalIntegrator& Integrate, mat& b, int vectorLevel=0)
    {
        SetLocalIntegrator(Integrate);
        cout<<"Integrating over "<<u_Internal.PhysicalGrpName<<"\n";
        NodePositions_u=std::vector<umat>(u_Internal.NoOfElementTypes);
        NodePositions_v=std::vector<umat>(u_Internal.NoOfElementTypes);
        for (int ElmntTyp=0; ElmntTyp<u_Internal.NoOfElementTypes; ElmntTyp++)
        {
            GetNodePostions(NodePositions_u[ElmntTyp], u_Internal.ElmntNodes[ElmntTyp], u_Internal.originalVctrLvl);
            GetNodePostions(NodePositions_v[ElmntTyp], v_Internal.ElmntNodes[ElmntTyp], v_Internal.originalVctrLvl);
            //cout<<"u_Internal.ElmntNodes[ElmntTyp] "<<u_Internal.ElmntNodes[ElmntTyp];
            //cout<<"NodePositions_u[ElmntTyp] "<<NodePositions_u[ElmntTyp]<<"\n";
            std::vector<mat> bTemp(numOfThreads);
            int ElmntDivision=u_Internal.NoOfElements[ElmntTyp]/numOfThreads;
            if(integrate->IsUsingFormMultiThread())
            {
                bTemp=std::vector<mat>(numOfThreads);
                std::thread Fill_A_Thread[numOfThreads];
                for (int thread=0; thread<numOfThreads; thread++)
                {
                    if(thread<numOfThreads-1)
                    {
                        Fill_A_Thread[thread]=std::thread(&SystemAssembler::IntegrateOverElements4Vector, this, std::ref(b),
                                                          std::ref(bTemp[thread]), thread, ElmntTyp,
                                                          ElmntDivision*thread, ElmntDivision*(thread+1));
                    }
                    else
                    {
                        Fill_A_Thread[thread]=std::thread(&SystemAssembler::IntegrateOverElements4Vector, this, std::ref(b),
                                                          std::ref(bTemp[thread]), thread, ElmntTyp,
                                                          ElmntDivision*thread, u_Internal.NoOfElements[ElmntTyp]);
                    }
                }
                for (int thread=0; thread<numOfThreads; thread++)
                {
                    Fill_A_Thread[thread].join();
                    b=b+bTemp[thread];
                }

            }
            else
            {
                bTemp=std::vector<mat>(1);
                IntegrateElements4Vector(b, bTemp[0], ElmntTyp);
                b=b+bTemp[0];
            }
        }
    }

    void RunScalarIntegration(GenericLocalIntegrator& Integrate, vec& I)
    {
        SetLocalIntegrator(Integrate);
        SetScalarSize(I);
        cout<<"Integrating over "<<u_Internal.PhysicalGrpName<<"\n";
        for (int ElmntTyp=0; ElmntTyp<u_Internal.NoOfElementTypes; ElmntTyp++)
        {
            int ElmntDivision=u_Internal.NoOfElements[ElmntTyp]/numOfThreads;
            if(integrate->IsUsingFormMultiThread())
            {
                //std::vector<double> I_Temp(numOfThreads);
                std::thread Fill_I_Thread[numOfThreads];
                for (int thread=0; thread<numOfThreads; thread++)
                {
                    if(thread<numOfThreads-1)
                    {
                        Fill_I_Thread[thread]=std::thread(&SystemAssembler::IntegrateOverElements4Scalar, this, std::ref(I), thread, ElmntTyp,
                                                          ElmntDivision*thread, ElmntDivision*(thread+1));
                    }
                    else
                    {
                        Fill_I_Thread[thread]=std::thread(&SystemAssembler::IntegrateOverElements4Scalar, this, std::ref(I), thread, ElmntTyp,
                                                          ElmntDivision*thread, u_Internal.NoOfElements[ElmntTyp]);
                    }
                }
                for (int thread=0; thread<numOfThreads; thread++)
                {
                    Fill_I_Thread[thread].join();
                }

            }
            else
            {
                for (int ElmntNmbr=0; ElmntNmbr<u_Internal.NoOfElements[ElmntTyp]; ElmntNmbr++)
                {
                    RunLocalIntegrationScalar();
                    I(ElmntNmbr)=integrate->GetResultingScalar();
                    integrate->GoToNextElement();
                }
             integrate->GoToNextElementType();
            }
        }
    }

private:
    FormMultiThread<GenericTrialFunction> a_Internal;
    GenericTrialFunction& u_Internal;
    TestFunctionGalerkin<GenericTrialFunction>& v_Internal;
    LocalIntegrator<GenericTrialFunction>* integrate;
    std::vector<umat> NodePositions_u;
    std::vector<umat> NodePositions_v;
    sp_mat* A_Internal;
    mat * b_Internal;
    int numOfThreads;

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

    void RunLocalIntegration(int thread)
    {
        integrate->local_intergrator(thread);
    }

    void RunLocalIntegrationVector()
    {
        integrate->local_intergrator_vector();
    }

    void RunLocalIntegrationVector(int thread)
    {
        integrate->local_intergrator_vector(thread);
    }

    void RunLocalIntegrationScalar()
    {
        integrate->local_integration_scalar();
    }

    void RunLocalIntegrationScalar(int thread)
    {
        integrate->local_integration_scalar(thread);
    }


    ///Responsible for adding over each element and then moving on to the next element for matrices. This is multithreaded.
    void IntegrateOverElements4Matrix(sp_mat &A, sp_mat &Atemp, int thread, int ElmntTyp, int ElmntStart, int ElmntEnd)
    {
        cout<<"Thread "<<thread<<" has been launched!!!\n";
        integrate->SetElementStartTo(thread, ElmntStart);
        //int thread=0;
        Atemp.set_size(A.n_rows, A.n_cols);
        for (int ElmntNmbr=ElmntStart; ElmntNmbr<ElmntEnd ; ElmntNmbr++)
        {
            RunLocalIntegration(thread);
            umat positions_u=NodePositions_u[ElmntTyp].row(integrate->GetElementNumber(thread));
            umat positions_v=NodePositions_v[ElmntTyp].row(integrate->GetElementNumber(thread));
            //cout<<"ElmntNmbr ="<<ElmntNmbr<<" ";
            sp_mat A_temp2=BatchFill_Atemp(A, ElmntTyp, positions_u, positions_v, thread);
            Atemp=Atemp+A_temp2;
            //a_Internal.NextElementNumber();
            integrate->GoToNextElement(thread);
        }
        cout<<"\n";
        //a_Internal.NextElementType();
        integrate->GoToNextElementType(thread);
    }

    ///Responsible for adding over each element and then moving on to the next element for matrices. This is single threaded.
    void IntegrateElement4Matrix(sp_mat& A, sp_mat& Atemp, int ElmntTyp)
    {
        Atemp.set_size(A.n_rows, A.n_cols);
        for (int ElmntNmbr=0; ElmntNmbr<u_Internal.NoOfElements[ElmntTyp] ; ElmntNmbr++)
        {
            RunLocalIntegration();
            umat positions_u=NodePositions_u[ElmntTyp].row(integrate->GetElementNumber());
            umat positions_v=NodePositions_v[ElmntTyp].row(integrate->GetElementNumber());
            //cout<<"ElmntNmbr ="<<ElmntNmbr<<"\n";
            sp_mat A_temp2=BatchFill_Atemp(A, ElmntTyp, positions_u, positions_v);
            Atemp=Atemp+A_temp2;
            //a_Internal.NextElementNumber();
            integrate->GoToNextElement();
        }
        //a_Internal.NextElementType();
        integrate->GoToNextElementType();
    }


    /// This function fill ups the A_temp or rather the local matrix for the element.
    /// The current value of the local element is stored in the 'a_Internal.ResultingMat' matrix.
    inline sp_mat BatchFill_Atemp(sp_mat& A, int& ElmntTyp, const umat& positions_u, const umat& positions_v, int thread=0)
    {
        //Configuration for Batch addition of matrix to global matrix
        bool add_values=false;
        bool sort_locations=true;
        bool check_for_zeros=true;
        //------------------------------------------------------------
        umat locations(2,positions_v.n_cols*positions_u.n_cols);
        //cout<<"positions_v.n_cols = "<<positions_v.n_cols<<" positions_u.n_cols = "<<positions_u.n_cols<<"\n";
        int locationPtr=0;
        for (int col=0;col<positions_u.n_cols;col++)
        {
            for (int row=0;row<positions_v.n_cols;row++)
            {
                umat locations1;
                locations1<<positions_v(row)<<endr<<positions_u(col)<<endr;
                locations.col(locationPtr)=locations1;
                locationPtr++;
            }
        }
        vec values=vectorise(mat(integrate->GetResultingMat(thread)));
        //cout<<mat(a_Internal.ResultingMat.t())<<"\n";
        //cout<<"Size of locations ="<<locations.n_rows<<", "<<locations.n_cols<<"\n";
        //cout<<"Size of values ="<<values.n_rows<<", "<<values.n_cols<<"\n";
        sp_mat A_temp=sp_mat(add_values, locations, values, A.n_rows, A.n_cols, sort_locations, check_for_zeros);
        return A_temp;
    }


    /// Responsible for adding over each element and then moving on to the next element for vectors. This is multithreaded.
    void IntegrateOverElements4Vector(mat &b, mat &btemp, int thread, int ElmntTyp, int ElmntStart, int ElmntEnd)
    {
        cout<<"Thread "<<thread<<" has been launched!!!\n";
        integrate->SetElementStartTo(thread, ElmntStart);
        btemp=zeros(b.n_rows, b.n_cols);
        //int thread=0;
        for (int ElmntNmbr=ElmntStart; ElmntNmbr<ElmntEnd; ElmntNmbr++)
        {
            RunLocalIntegrationVector(thread);
            umat positions_v=NodePositions_v[ElmntTyp].row(integrate->GetElementNumber(thread));
            //cout<<"Neumann Condition applied over "<<u_Internal.Msh->NodalCoordinates.rows(temp)<<"\n";
            umat positions2={0};
            btemp.submat(positions_v, positions2)=btemp.submat(positions_v, positions2)+integrate->GetResultingVector(thread);
            //a_Internal.NextElementNumber();
            integrate->GoToNextElement(thread);
        }
        //a_Internal.NextElementType();
        integrate->GoToNextElementType(thread);
    }

    /// Responsible for adding over each element and then moving on to the next element for vectors. This is single threaded.
    void IntegrateElements4Vector(mat &b, mat &btemp, int ElmntTyp)
    {
        btemp=zeros(b.n_rows, b.n_cols);
        for (int ElmntNmbr=0; ElmntNmbr<u_Internal.NoOfElements[ElmntTyp]; ElmntNmbr++)
        {
            RunLocalIntegrationVector();
            umat positions_v=NodePositions_v[ElmntTyp].row(integrate->GetElementNumber());
            //cout<<"Neumann Condition applied over "<<u_Internal.Msh->NodalCoordinates.rows(temp)<<"\n";
            umat positions2={0};
            btemp.submat(positions_v, positions2)=btemp.submat(positions_v, positions2)+integrate->GetResultingVector();
            //a_Internal.NextElementNumber();
            integrate->GoToNextElement();
        }
        //a_Internal.NextElementType();
        integrate->GoToNextElementType();
    }

    void IntegrateOverElements4Scalar(vec &I, int thread, int ElmntTyp, int ElmntStart, int ElmntEnd)
    {
        integrate->SetElementStartTo(thread, ElmntStart);
        //cout<<"Size of I = ("<<I.n_rows<<", "<<I.n_cols<<")\n";
        for (int ElmntNmbr=ElmntStart; ElmntNmbr<ElmntEnd; ElmntNmbr++)
        {
            RunLocalIntegrationScalar(thread);
            I(ElmntNmbr,0)=integrate->GetResultingScalar(thread);
            integrate->GoToNextElement(thread);
        }
        integrate->GoToNextElementType(thread);
        cout<<"\n";
    }

};

#endif // SYSTEMASSEMBLY_HPP
