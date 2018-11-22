 function [pos NodalCoordinates] = GetNodeData(FileContents, pos , WordStart, WordEnd)
    pos=pos+1;
    [strng]=GetWord(pos, FileContents, WordStart, WordEnd);
    pos=pos+1;
    NumberOfNodes=str2num(strng);
    for i=1:NumberOfNodes
        [strng]=GetWord(pos, FileContents, WordStart, WordEnd);
        pos=pos+1;
        NodeNumber=str2num(strng);
        for j=1:3
            [strng]=GetWord(pos, FileContents, WordStart, WordEnd);
            pos=pos+1;
            NodalCoordinates(i,j)=str2num(strng);
        end
    end
end
