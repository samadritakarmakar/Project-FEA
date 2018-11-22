%function MeshReader(FileName)
clear
FileName=input('Enter the Mesh File to read: ','s');
    FileContents=fileread(FileName);
    FileLength=length(FileContents);
    SpacePosition=findstr(FileContents,' ');
    NewLinePosition=findstr(FileContents,sprintf('\n'));
    InBtwnLetterPos=[SpacePosition, NewLinePosition];
    InBtwnLetterPos=sort(InBtwnLetterPos);
    if(InBtwnLetterPos(end)!=FileLength)
       InBtwnLetterPos(end)=[InBtwnLetterPos FileLength+1]; 
    end
    WordStart=InBtwnLetterPos(1:end-1)+ones(size(InBtwnLetterPos(1:end-1)));
    WordStart=[1 WordStart];
    WordEnd=InBtwnLetterPos(1:end)-ones(size(InBtwnLetterPos(1:end)));
    for i=1:length(WordEnd)
        Word=GetWord(i, FileContents, WordStart, WordEnd);
        switch(Word)
        case '$MeshFormat'
            [i Version]=GetVersion(i, FileContents, WordStart, WordEnd);
        case '$Nodes'
            [i, NodalCoordinates]=GetNodeData(FileContents, i, WordStart, WordEnd);
        case '$Elements'
            [i, ElementData]=GetElementData(FileContents, i, WordStart, WordEnd);
        end
    end
%end
