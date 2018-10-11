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
Property.type=ElementData.Type;
Property.degree=ElementData.Degree;
[NumOfElements , ~]=size(ElementData.ElementNodes);

file = LoadGaussFile(Property);

load (file)
[GaussLength ~]=size(data);
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
NodePositons=GetNodePostions(ElementData.ElementNodes,vectorLevel);

for ElementNum=1:NumOfElements
    NodalCoord=NodalCoordinates (ElementData.ElementNodes(ElementNum,:),:);
    xCoords=NodalCoord(:,1);
    yCoords=NodalCoord(:,2);
    zCoords=NodalCoord(:,3);
    
    OtherData.Property=Property;
    OtherData.xCoords=xCoords;
    OtherData.yCoords=yCoords;
    OtherData.zCoords=zCoords;

    if (strcmp(Property.type,'1D'))
        for GaussPoint=1:GaussLength
            epsilon=data(GaussPoint,3:end);
            w=data(GaussPoint,2);
            phi=ShapeFunction(epsilon, Property);
            %have doubts when using higher order elements
            xDiff=xCoords(end)-xCoords(1);
            yDiff=yCoords(end)-yCoords(1);
            zDiff=zCoords(end)-zCoords(1);
            if (vectorLevel==1)
                ElementLength=(xDiff);
                cos_xyz=1;
            elseif (vectorLevel==2)
                ElementLength=sqrt(xDiff^2+yDiff^2);
                cos_xyz=[xDiff,yDiff]*(1/ElementLength);
            elseif (vectorLevel==3)
                ElementLength=sqrt(xDiff^2+yDiff^2+zDiff^2);
                cos_xyz=[xDiff,yDiff,zDiff]*(1/ElementLength);
            end
            cosMatrix=zeros(ElementData.NumOfElementNodes,vectorLevel*ElementData.NumOfElementNodes);
            xVector=zeros(vectorLevel*ElementData.NumOfElementNodes,1);
            for Row=1:ElementData.NumOfElementNodes;
                if (vectorLevel==1)
                   Coords= [xCoords(Row)];
                elseif (vectorLevel==2)
                   Coords= [xCoords(Row);yCoords(Row)];
                elseif (vectorLevel==3)
                    Coords= [xCoords(Row);yCoords(Row);zCoords(Row)];
                end
                cosMatrix(Row,(vectorLevel*Row-(vectorLevel-1)):(vectorLevel*Row))=cos_xyz;  %A First Course in the Finite Element Method - Daryl L. Logan eqn 3.7.8
                xVector((vectorLevel*Row-(vectorLevel-1)):(vectorLevel*Row),1)= Coords;
            end
            OtherData.xCoords=cosMatrix*xVector;
            x=interpolateX(epsilon, OtherData);
            u=phi;
            v=phi';
            F=jacobian(@interpolateX, epsilon, OtherData); %Jonathan Whiteley Finite Element Methods A Practical Guide eqn 7.39
            duBYdx=jacobian(@ShapeFunction,epsilon,Property)*inv(F); %Jonathan Whiteley Finite Element Methods A Practical Guide eqn 7.39
            duBYdx=duBYdx';
            dvBydx=duBYdx;
            [LHSmatrixLocalGauss RHSmatrixLocalGauss RHSvectorLocalGauss]=UserFunction(x,u,v,duBYdx,dvBydx, ElementNum);
            if(isAvector)
                LHSmatrixLocalGauss=cosMatrix'*LHSmatrixLocalGauss*cosMatrix; %A First Course in the Finite Element Method - Daryl L. Logan eqn 3.7.8
                RHSvectorLocalGauss=cosMatrix'*RHSvectorLocalGauss; %T'*f Change local matrix to global matrix; A First Course in the Finite Element Method  3.4.16
            end
            if (GaussPoint==1)
                LHSmatrixLocal=zeros(size(LHSmatrixLocalGauss));
                RHSmatrixLocal=zeros(size(RHSmatrixLocalGauss));
                RHSvectorLocal=zeros(size(RHSvectorLocalGauss));
            end
             %Jonathan Whiteley Finite Element Methods A Practical Guide eqn 7.36
            LHSmatrixLocal=w*LHSmatrixLocalGauss*det(F) + LHSmatrixLocal;
            RHSmatrixLocal=w*RHSmatrixLocalGauss*det(F) + RHSmatrixLocal;
            RHSvectorLocal=w*RHSvectorLocalGauss*det(F) + RHSvectorLocal;
        end
        if (ElementNum==1)
            LHSmatrix=zeros(dof);
            RHSvector=zeros(dof,1);
            if (length(RHSmatrixLocalGauss)==0)
                RHSmatrix=[];
                1
            else
                RHSmatrix=zeros(dof);
            end
        end
    elseif(strcmp(Property.type, 'Triangle'))
        for GaussPoint=1:GaussLength
            epsilon=data(GaussPoint,3:end);
            phi = ShapeFunction(epsilon, Property);
            x=interpolateX(epsilon, OtherData);
            y=interpolateY(epsilon, OtherData);
            u=phi;
            v=phi';
        end
    elseif(strcmp(Property.type, 'Quadrilateral'))
        for GaussPoint1=1:GaussLength
            for GaussPoint2=1:GaussLength
                epsilon=[data(GaussPoint1,3:end),data(GaussPoint2,3:end)];
                phi = ShapeFunction(epsilon, Property);
                x=interpolateX(epsilon, OtherData);
                y=interpolateY(epsilon, OtherData);
                u=phi;
                v=phi';
            end
        end
    end
    LHSmatrix(NodePositons(ElementNum,:),NodePositons(ElementNum,:))=LHSmatrix(NodePositons(ElementNum,:),NodePositons(ElementNum,:))+ LHSmatrixLocal;
    if (length(RHSmatrix)~=0)
        RHSmatrix(NodePositons(ElementNum,:),NodePositons(ElementNum,:))=RHSmatrix(NodePositons(ElementNum,:),NodePositons(ElementNum,:))+ RHSmatrixLocal;
    end
    RHSvector(NodePositons(ElementNum,:),1)=RHSvector(NodePositons(ElementNum,:),1)+RHSvectorLocal;
end


