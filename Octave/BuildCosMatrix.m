function [cosMatrix xVector]= BuildCosMatrix(xCoords, yCoords, zCoords, ElementData, vectorLevel)
            xDiff=xCoords(end)-xCoords(1);
            yDiff=yCoords(end)-yCoords(1);
            zDiff=zCoords(end)-zCoords(1);
            DiffVector=[xDiff;yDiff;zDiff];
           %Old implementation retained for understading new implementation below. 
           % if (vectorLevel==1)
           %     ElementLength=(xDiff);
           %     cos_xyz=1;
           % elseif (vectorLevel==2)
           %     ElementLength=sqrt(xDiff^2+yDiff^2);
           %     cos_xyz=[xDiff,yDiff]*(1/ElementLength);
           % elseif (vectorLevel==3)
           %     ElementLength=sqrt(xDiff^2+yDiff^2+zDiff^2);
           %     cos_xyz=[xDiff,yDiff,zDiff]*(1/ElementLength);
           % end

            SquareOfDiffVector=DiffVector.^2;
            SumOfSquare=sum(SquareOfDiffVector(1:vectorLevel),1); %vector Squared only upto the specified Vector Level
            ElementLength=sqrt(SumOfSquare);
            cos_xyz=DiffVector(1:vectorLevel)'./ElementLength; %vector Saved only upto the specified Vector Level
            cosMatrix=zeros(ElementData.NumOfElementNodes,vectorLevel*ElementData.NumOfElementNodes);
            xVector=zeros(vectorLevel*ElementData.NumOfElementNodes,1);
            CoordsTempStore=[xCoords, yCoords, zCoords]; %Temporary storage of coordiantes.
            for Row=1:ElementData.NumOfElementNodes;
                Coords=CoordsTempStore(Row,1:vectorLevel)'; %vector Saved only upto the specified Vector Level
                cosMatrix(Row,(vectorLevel*Row-(vectorLevel-1)):(vectorLevel*Row))=cos_xyz;  %A First Course in the Finite Element Method - Daryl L. Logan eqn 3.7.8; [... 0 0 Cx Cy Cz 0 0 ...]
                xVector((vectorLevel*Row-(vectorLevel-1)):(vectorLevel*Row),1)= Coords;
            end
end
