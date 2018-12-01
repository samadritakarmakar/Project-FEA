function [pos ElementData] = GetElementData(FileContents, pos, WordStart, WordEnd)
    GmshElementNumBuffer=[];
    pos=pos+1
    [strng]=GetWord(pos, FileContents, WordStart, WordEnd);
    loopEnd=str2num(strng);
    j=0;
    for i=1:loopEnd
        [strng]=GetWord(pos, FileContents, WordStart, WordEnd);
        pos=pos+1;
        [strng]=GetWord(pos, FileContents, WordStart, WordEnd);
        pos=pos+1;
        GmshElementNum=str2num(strng)
        [Type Degree NumOfElementNodes NumofDimensions Supported]=GetElemProperty(GmshElementNum);
        if (Supported)
            if (GmshElementNum~=GmshElementNumBuffer)
                j=j+1;
                l=0;
            end
            l=l+1;
            ElementType(j).Type = Type;
            ElementType(j).Degree = Degree;
            ElementType(j).NumOfElementNodes = NumOfElementNodes;
            ElementType(j).NumofDimensions = NumofDimensions;
            [strng pos]=SearchNextString(FileContents, pos, FileLength);
            GmshNumOfTags=str2num(strng);
            pos = 2*GmshNumOfTags+pos;
            for k=1:NumOfElementNodes
                [strng]=GetWord(pos, FileContents, WordStart, WordEnd);
                pos=pos+1;
                ElementType(j).ElementNodes(l,k)=str2num(strng);
            end
        else
          %  pos = SkipToNextLine(FileContents, pos, FileLength);
        end
    GmshElementNumBuffer=GmshElementNum;
    end
    fprintf('The supported Elements found are\n');
    for l=1:j
        fprintf('%i %s \n',l, ElementType(l).Type);
    end
    n=input('Enter the Element to integrate over: ');
    ElementData=ElementType(n);
    ElementData.ContainsNodes=SortAndKeep(ElementData.ElementNodes);
end
        
        
