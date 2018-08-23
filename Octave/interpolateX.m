function [x]=interpolateX(epsilon, OtherData)
    Property=OtherData.Property;
    xCoords=OtherData.xCoords;
    yCoords=OtherData.yCoords;
    zCoords=OtherData.zCoords;
    phi=ShapeFunction(epsilon, Property);
    x=phi*xCoords;
    y=phi*yCoords;
    z=phi*zCoords;
end
