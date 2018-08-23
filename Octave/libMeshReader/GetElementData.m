function ElementData = GetElementData(FileContents, pos, FileLength)
    [strng pos]=SearchNextString(FileContents, pos, FileLength);
    loopEnd=str2num(strng);
    j=0;
    for i=1:loopEnd
        [strng pos]=SearchNextString(FileContents, pos, FileLength);
        [strng pos]=SearchNextString(FileContents, pos, FileLength);
        GmshElementNum=str2num(strng);
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
                [strng pos]=SearchNextString(FileContents, pos, FileLength);
                ElementType(j).ElementNodes(l,k)=str2num(strng);
            end
        else
            pos = SkipToNextLine(FileContents, pos, FileLength);
        end
    GmshElementNumBuffer=GmshElementNum;
    end
    fprintf('The supported Elements found are\n');
    for l=1:j
        fprintf('%i %s \n',l, ElementType(l).Type);
    end
    n=input('Enter the Element to integrate over: ');
    ElementData=ElementType(n);
end       
        

function [Type Degree NumOfElementNodes NumofDimensions Supported]=GetElemProperty(GmshElementNum)
    switch(GmshElementNum)
        case  1
            Type='1D';
            Degree='1';
            NumOfElementNodes=2;
            NumofDimensions=1;
            Supported=1;
        case 2
            Type='Triangle';
            Degree='1';
            NumOfElementNodes=3;
            NumofDimensions=2;
            Supported=1;
        case 3
            Type='Quadrilateral';
            Degree='1';
            NumOfElementNodes=4;
            NumofDimensions=2;
            Supported=1;
        case 4
            Type='Tetrahedral';
            Degree='1';
            NumOfElementNodes=4;
            NumofDimensions=3;
            Supported=1;
        case 5
            Type='Hexahedral'
            Degree='1';
            NumOfElementNodes=8;
            NumofDimensions=3;
            Supported=1;
        case 8
            Type='1D';
            Degree='2';
            NumOfElementNodes=3;
            NumofDimensions=1;
            Supported=1;
        case 9
            Type='Triangle';
            Degree='2';
            NumOfElementNodes=6;
            NumofDimensions=2;
            Supported=1;
        case 10
            Type='Quadrilateral';
            Degree='Biquadratic';
            NumOfElementNodes=9;
            NumofDimensions=2;
            Supported=1;
        case 11
            Type='Tetrahedral';
            Degree='2';
            NumOfElementNodes=10;
            NumofDimensions=3;
            Supported=1;
        case 16
            Type='Quadrilateral';
            Degree='Serendipity';
            NumOfElementNodes=8;
            NumofDimensions=2;
            Supported=1;
        case 21
            Type='Triangle';
            Degree='3';
            NumOfElementNodes=10;
            NumofDimensions=2;
            Supported=1;
        case 26
            Type='1D';
            Degree='3';
            NumOfElementNodes=4;
            NumofDimensions=1;
            Supported=1;
        otherwise
            Type='';
            Degree='';
            NumOfElementNodes=0;
            NumofDimensions=0;
            Supported=0;
    end
end
    
function pos=SkipToNextLine(FileContents, pos, FileLength)
    do
        strng = FileContents(pos);
        if (~strcmp(strng,sprintf('\n')))
            pos=pos+1;
        end
    until(FileContents(pos)==sprintf('\n') || pos==FileLength)
end
