clear;
fprintf('A general FEA solver with .msh file and weak form as basic inputs\n');
fprintf('    Copyright (C) 2018  Samadrita Karmkar (samadritakarmakar@gmail.com)\n\n ');
 
fprintf('    This program is free software: you can redistribute it and/or modify\n ');
fprintf('    it under the terms of the GNU General Public License as published by\n ');
fprintf('    the Free Software Foundation, either version 3 of the License, or\n ');
fprintf('    (at your option) any later version.\n\n ');

fprintf('    This program is distributed in the hope that it will be useful,\n ');
fprintf('    but WITHOUT ANY WARRANTY; without even the implied warranty of\n ');
fprintf('    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n ');
fprintf('    GNU General Public License for more details.\n\n')

fprintf('    You should have received a copy of the GNU General Public License\n \n');
fprintf('    along with this program.  If not, see <https://www.gnu.org/licenses/>.\n\n\n');
    


addpath('libMeshReader/')
FileName=input('Enter file name of Mesh file: ','s');
[ElementData NodalCoordinates] = MeshReader(FileName);
Property.Type=ElementData.Type;
Property.degree=ElementData.Degree;
NumOfElements=size(ElementData.ElementNodes,1);
%loads the Gauss File of the correct Poperty

file = LoadGaussFile(Property);
load (file, '-mat');
[GaussLength DependentEpsilon]=size(data);
DependentEpsilon=DependentEpsilon-2; %Factoring for the numbering and the weights in the Gauss Files
isAvectorCh=input('Is the quantity in question a vector? (y/n) ','s');
if (strcmp(isAvectorCh,'y')||strcmp(isAvectorCh,'Y'))
    isAvector=1;
    vectorLevel=input('Enter vector Level. 2 for 2D and 3 for 3D: ');
    dof=vectorLevel*length(ElementData.ContainsNodes); %Here it the total number of nodes 
else
    isAvector=0;
    vectorLevel=1;
    dof=1*length(ElementData.ContainsNodes);
end
%Get Node Positions in the Matrices and the Vectors
NodePositons=GetNodePostions(ElementData.ElementNodes,vectorLevel); %Positions in the matrix and vectors which are to be filled.
if (strcmp(Property.Type,'1D'))
GaussLength1=GaussLength;
GaussLength2=1;
GaussLength3=1;
NoOfInterpolatedCoords=1;
NoOfIndependentEpsilon=1;
elseif(strcmp(Property.Type, 'Triangle'))
GaussLength1=GaussLength;
GaussLength2=1;
GaussLength3=1;
NoOfInterpolatedCoords=2;
NoOfIndependentEpsilon=1;
elseif(strcmp(Property.Type, 'Quadrilateral'))
GaussLength1=GaussLength;
GaussLength2=GaussLength;
GaussLength3=1;
NoOfInterpolatedCoords=2;
NoOfIndependentEpsilon=2;
elseif(strcmp(Property.Type, 'Tetrahedral'))
GaussLength1=GaussLength;
GaussLength2=1;
GaussLength3=1;
NoOfInterpolatedCoords=3;
NoOfIndependentEpsilon=1;
elseif(strcmp(Property.Type, 'Hexahedral'))
GaussLength1=GaussLength;
GaussLength2=GaussLength;
GaussLength3=GaussLength;
NoOfInterpolatedCoords=3;
NoOfIndependentEpsilon=3;
end

%Added to test and improve code execution time
profile on
tic    

for ElementNum=1:NumOfElements
    NodalCoord=NodalCoordinates (ElementData.ElementNodes(ElementNum,:),:);
    xCoords=NodalCoord(:,1);
    yCoords=NodalCoord(:,2);
    zCoords=NodalCoord(:,3);
    
    OtherData.Property=Property;
    OtherData.xCoords=xCoords;
    OtherData.yCoords=yCoords;
    OtherData.zCoords=zCoords;
    OtherData.NoOfInterpolatedCoords=NoOfInterpolatedCoords;
%When considering only 1D elements
  %  if (strcmp(Property.Type,'1D'))
    for GaussPoint1=1:GaussLength1
        for GaussPoint2=1:GaussLength2
            for GaussPoint3=1:GaussLength3
                epsilonTemp=[data(GaussPoint1,3:end),data(GaussPoint2,3:end),data(GaussPoint3,3:end)];
                epsilon=epsilonTemp(1:NoOfIndependentEpsilon*DependentEpsilon);
                wTemp=[data(GaussPoint1,2),data(GaussPoint2,2),data(GaussPoint3,2)];
                w=[wTemp(1:NoOfIndependentEpsilon),ones(1,3-NoOfIndependentEpsilon)];
                switch Property.Type
                case '1D'
                    %Doubtful about Higher Order Elements.
                    %Done for Transformation of Vectors.
                    [cosMatrix xVector]= BuildCosMatrix(xCoords, yCoords, zCoords, ElementData, vectorLevel);
                    %A First Course in the Finite Element Method - Daryl L. Logan eqn 3.7.6
                    OtherData.xCoords=cosMatrix*xVector; %Translate vector [x] containing x y and z coordinates 
                end
                phiTemp=ShapeFunction(epsilon, Property);%[phi1,phi2,phi3,.....]
                xInterpolated=interpolateX(epsilon, OtherData);
                VectorizedX=VectorizeX(OtherData, vectorLevel);
                %u=VectorizePhi(phiTemp, NoOfInterpolatedCoords);
                u=VectorizePhi(phiTemp, vectorLevel);
                v=u';
                F=jacobian(@interpolateX, epsilon, OtherData); %Jonathan Whiteley Finite Element Methods A Practical Guide eqn 7.39
                graduTemp=jacobian(@ShapeFunction,epsilon,Property)*inv(F); %Jonathan Whiteley Finite Element Methods A Practical Guide eqn 7.39
                %gradu=VectorizeGradu(graduTemp, NoOfInterpolatedCoords);%[phi1/dx phi1/dy  phi1/dz; phi2/dx  phi2/dy  phi2/z; phi3/dx phi3/dy phi3/dz;  ...]
                gradu=VectorizeGradu(graduTemp, vectorLevel);% Finite Element Methods for Fluid Problems- Jean Donea and Antonio Huerta eqn 6.24
                gradv=gradu;
                gradu=gradu';
                [LHSmatrixLocalGauss RHSmatrixLocalGauss RHSvectorLocalGauss]=UserFunction(xInterpolated,VectorizedX,u,v,gradu,gradv, ElementNum);
                switch (isAvector && strcmp(Property.Type,'1D'))
                    case 1
                        LHSmatrixLocalGauss=cosMatrix'*LHSmatrixLocalGauss*cosMatrix; %A First Course in the Finite Element Method - Daryl L. Logan eqn 3.7.8
                        RHSvectorLocalGauss=cosMatrix'*RHSvectorLocalGauss; %T'*f Change local matrix to global matrix; A First Course in the Finite Element Method  3.4.16
                end
                LocalGaussMatrices.LHSmatrixLocalGauss=LHSmatrixLocalGauss;
                LocalGaussMatrices.RHSmatrixLocalGauss=RHSmatrixLocalGauss;
                LocalGaussMatrices.RHSvectorLocalGauss=RHSvectorLocalGauss;

                switch GaussPoint1
                case 1
                    LHSmatrixLocal=zeros(size(LHSmatrixLocalGauss));
                    RHSmatrixLocal=zeros(size(RHSmatrixLocalGauss));
                    RHSvectorLocal=zeros(size(RHSvectorLocalGauss));
                end
            
                LocalMatrices.LHSmatrixLocal=LHSmatrixLocal;
                LocalMatrices.RHSmatrixLocal=RHSmatrixLocal;
                LocalMatrices.RHSvectorLocal=RHSvectorLocal;
                [LHSmatrixLocal RHSmatrixLocal RHSvectorLocal] = BuildLocalMatrices(LocalGaussMatrices, LocalMatrices, F, w);
            end
        end
    end
    switch ElementNum
    case 1
        LHSmatrix=spalloc(dof,dof,NumOfElements*ElementData.NumOfElementNodes*vectorLevel);
        %LHSmatrix=zeros(dof);
        %RHSvector=zeros(dof,1);
        switch length(RHSmatrixLocalGauss)
        case 0
            RHSmatrix=[];
        otherwise
            %RHSmatrix=zeros(dof);
            RHSmatrix=spalloc(dof,dof,NumOfElements*ElementData.NumOfElementNodes*vectorLevel);
        end
        switch length(RHSvectorLocalGauss)
        case 0
            RHSvector=[];
        otherwise
            RHSvector=zeros(dof,1);
        end
    end
    LHSmatrix(NodePositons(ElementNum,:),NodePositons(ElementNum,:))=LHSmatrix(NodePositons(ElementNum,:),NodePositons(ElementNum,:))+ LHSmatrixLocal;
    if (length(RHSmatrix)~=0)
        RHSmatrix(NodePositons(ElementNum,:),NodePositons(ElementNum,:))=RHSmatrix(NodePositons(ElementNum,:),NodePositons(ElementNum,:))+ RHSmatrixLocal;
    end
    RHSvector(NodePositons(ElementNum,:),1)=RHSvector(NodePositons(ElementNum,:),1)+RHSvectorLocal;
end
timeofExection=toc
profile off
