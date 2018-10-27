function gradu=VectorizeGradu(graduTemp, NoOfInterpolatedCoords)
    NoOfGraduTempRows=size(graduTemp,1);
    NoOfGraduTempColumns=size(graduTemp,2);
    gradu=spalloc(NoOfInterpolatedCoords*NoOfGraduTempRows, NoOfInterpolatedCoords*NoOfGraduTempColumns, NoOfGraduTempRows*NoOfGraduTempColumns);
    j=1;
    k=1;
    for i=1:NoOfInterpolatedCoords*NoOfGraduTempRows;
        startingColumn=k*NoOfInterpolatedCoords-(NoOfInterpolatedCoords-1); 
        gradu(i,startingColumn:startingColumn+(NoOfGraduTempColumns-1))=graduTemp(j,:);
        j=j+(mod(i,NoOfInterpolatedCoords)==0)*1;%Increases by one only when 'i' is divisible by NoOfInterpolatedCoords;
        k=k*(mod(i,NoOfInterpolatedCoords)~=0);%Becomes zero when 'i' is divisible by NoOfInterpolatedCoords; Forces k to start from zero again.
        k=k+1;
    end
end
