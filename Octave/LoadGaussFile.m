function file = LoadGaussFile(Property)
    if (strcmp(Property.type,'1D'))
        GaussDegree=2*str2num(Property.degree);
        file=strcat('Data/n',num2str(GaussDegree));
        elseif (strcmp(Property.type,'Triangle'))
            switch(Property.degree)
                case '1'
                    GaussDegree=1;
                case '2'
                    GaussDegree=3;
                case '3'
                    GaussDegree=4;
            end
        file=strcat('DataTriangle/n',num2str(GaussDegree));
    elseif (strcmp(Property.type,'Quadrilateral'))
        switch(Property.degree)
            case '1'
                GaussDegree=2;
            case 'Biquadratic'
                GaussDegree=3;
            case 'Serendipity'
                GaussDegree=3;
        end
        file=strcat('Data/n',num2str(GaussDegree));
    elseif (strcmp(Property.type,'Tetrahedral'))
        switch(Property.degree)
            case '1'
                GaussDegree=1;
            case '2'
                GaussDegree=4;
            case '3'
                GaussDegree=5;
        end
        file=strcat('DataTetrahedral/n',num2str(GaussDegree))
    elseif (strcmp(Property.type,'Hexahedral'))
        switch(Property.degree)
            case '1'
                GaussDegree=2;
            case '2'
                GaussDegree=3;
            case '3'
                GaussDegree=5;
        end 
        file=strcat('Data/n',num2str(GaussDegree));
    end
end 
