function VectorizedX=VectorizeX(OtherData, vectorLevel)
    Property=OtherData.Property;
    xCoords=OtherData.xCoords;
    yCoords=OtherData.yCoords;
    zCoords=OtherData.zCoords;
    for i=1:length(xCoords)
        VectorizedX_Temp=[xCoords(i);yCoords(i);zCoords(i)];
        VectorizedX(vectorLevel*i-(vectorLevel-1):vectorLevel*i,1)=VectorizedX_Temp(1:vectorLevel,1);
    end
end
