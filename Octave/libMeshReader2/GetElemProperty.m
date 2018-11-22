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
    

