function file = LoadGaussFile(Property)
    if (strcmp(Property.Type,'1D'))
        GaussDegree=2*str2num(Property.degree);
        file=strcat('Data2/n',num2str(GaussDegree));
    elseif (strcmp(Property.Type,'Triangle'))
        switch(Property.degree)
            case '1'
                GaussDegree=1;
            case '2'
                GaussDegree=3;
            case '3'
                GaussDegree=4;
        end
        file=strcat('DataTriangle2/n',num2str(GaussDegree));
    elseif (strcmp(Property.Type,'Quadrilateral'))
        switch(Property.degree)
            case '1'
                GaussDegree=2;
            case 'Biquadratic'
                GaussDegree=3;
            case 'Serendipity'
                GaussDegree=3;
        end
        file=strcat('Data2/n',num2str(GaussDegree));
    elseif (strcmp(Property.Type,'Tetrahedral'))
        switch(Property.degree)
            case '1'
                GaussDegree=1;
            case '2'
                GaussDegree=4;
            case '3'
                GaussDegree=5;
        end
        file=strcat('DataTetrahedral2/n',num2str(GaussDegree));
    elseif (strcmp(Property.Type,'Hexahedral'))
        switch(Property.degree)
            case '1'
                GaussDegree=2;
            case '2'
                GaussDegree=3;
            case '3'
                GaussDegree=5;
        end 
        file=strcat('Data2/n',num2str(GaussDegree));
    end
end 
