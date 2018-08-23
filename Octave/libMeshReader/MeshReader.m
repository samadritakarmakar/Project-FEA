function [ElementData NodalCoordinates] = MeshReader(FileName)
    FileContents=fileread(FileName);
    FileLength=length(FileContents);
    pos=1;
    while (pos<FileLength)
        [strng pos]=SearchNextString(FileContents, pos, FileLength);
        pos=pos+1;
        if (strcmp(strng,'$MeshFormat'))
            [strng pos]=SearchNextString(FileContents, pos, FileLength);
            fprintf('The Gmsh Format is %s\n',strng);
            if (str2num(strng)>2.2 || str2num(strng)<2)
                fprintf('The format is unsupported\n');
            end
        elseif (strcmp(strng, '$Nodes'))
            NodalCoordinates = GetNodeData(FileContents, pos, FileLength);
        elseif (strcmp(strng, '$Elements'))
            ElementData = GetElementData(FileContents, pos, FileLength);
        end
    end
end

