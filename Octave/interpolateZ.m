function [z]=interpolateZ(epsilon, OtherData)
    Property=OtherData.Property;
    zCoords=OtherData.zCoords;
    phi=ShapeFunction(epsilon, Property);
    z=phi*zCoords;
end
