function [xyz]=interpolateX(epsilon, OtherData)
    Property=OtherData.Property;
    xCoords=OtherData.xCoords;
    yCoords=OtherData.yCoords;
    zCoords=OtherData.zCoords;
    NoOfInterpolatedCoords=OtherData.NoOfInterpolatedCoords;
    phi=ShapeFunction(epsilon, Property);
    x=phi*xCoords;
    y=phi*yCoords;
    z=phi*zCoords;
    xyz=zeros(1,NoOfInterpolatedCoords);
    xyz_temporary=[x;y;z];
    xyz=xyz_temporary(1:NoOfInterpolatedCoords);
end
