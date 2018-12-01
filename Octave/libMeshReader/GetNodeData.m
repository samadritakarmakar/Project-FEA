 function NodalCoordinates = GetNodeData(FileContents, pos, FileLength)
    [strng pos]=SearchNextString(FileContents, pos, FileLength);
    NumberOfNodes=str2num(strng);
    for i=1:NumberOfNodes
        [strng pos]=SearchNextString(FileContents, pos, FileLength);
        NodeNumber=str2num(strng);
        for i=1:3
            [strng pos]=SearchNextString(FileContents, pos, FileLength);
            NodalCoordinates(NodeNumber,i)=str2num(strng);
        end
    end
end
