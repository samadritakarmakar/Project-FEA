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
            phi=ShapeFunction(epsilon, Property);
            xDiff=xCoords(end)-xCoords(1);
            yDiff=yCoords(end)-yCoords(1);
            zDiff=zCoords(end)-zCoords(1);
            ElementLength=sqrt(xDiff^2+yDiff^2+zDiff^2);
            cos_xyz=[xDiff,yDiff,zDiff]*(1/ElementLength);
            cosMatrix=zeros(ElementData.NumOfElementNodes,3*ElementData.NumOfElementNodes);
            xVector=zeros(3*ElementData.NumOfElementNodes,1);
            for Row=1:ElementData.NumOfElementNodes;
                cosMatrix(Row,(3*Row-2):(3*Row))=cos_xyz;
                xVector((3*Row-2):(3*Row),1)= [xCoords(Row);yCoords(Row);zCoords(Row)];
            end
            OtherData.xCoords=cosMatrix*xVector;
            x=interpolateX(epsilon, OtherData);
            u=phi;
            v=phi';
            duBYdx=jacobian(@ShapeFunction,epsilon,Property)*inv(jacobian(@interpolateX, epsilon, OtherData));
            duBYdx=duBYdx';
            dvBydx=duBYdx;
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
end


