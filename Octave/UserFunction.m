function [LHSmatrix RHSmatrix RHSvector]=UserFunction(xInterpolated,VectorizedX,u,v,gradu,gradv, ElementNum);

LHSmatrix = InnerProduct(gradv,gradu);
RHSmatrix=[];
RHSvector=Dot(v,f(VectorizedX));
