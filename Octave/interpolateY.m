function [y]=interpolateY(epsilon, OtherData)
    Property=OtherData.Property;
    yCoords=OtherData.yCoords;
    phi=ShapeFunction(epsilon, Property);
    y=phi*yCoords;
end
