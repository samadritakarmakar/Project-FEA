function gradu=VectorizeGradu(graduTemp, vectorLevel)
    NoOfGraduTempRows=size(graduTemp,1);
    NoOfGraduTempColumns=size(graduTemp,2);
    gradu=spalloc(vectorLevel*NoOfGraduTempRows, vectorLevel*NoOfGraduTempColumns, NoOfGraduTempRows*NoOfGraduTempColumns);
    j=1;
    k=1;
    for i=1:vectorLevel*NoOfGraduTempRows;
        startingColumn=k*vectorLevel-(vectorLevel-1); 
        gradu(i,startingColumn:startingColumn+(NoOfGraduTempColumns-1))=graduTemp(j,:);
        j=j+(mod(i,vectorLevel)==0)*1;%Increases by one only when 'i' is divisible by vectorLevel;
        k=k*(mod(i,vectorLevel)~=0);%Becomes zero when 'i' is divisible by vectorLevel; Forces k to start from zero again.
        k=k+1;
    end
end
