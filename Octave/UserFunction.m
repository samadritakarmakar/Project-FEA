function [LHSmatrix RHSmatrix RHSvector]=UserFunction(xInterpolated,VectorizedX,u,v,gradu,gradv, ElementNum);

LHSmatrix = ScalarDot(gradv,gradu);
RHSmatrix=[];
RHSvector=ScalarDot(v,u*f(VectorizedX)); %u added to enable integration at Gauss Point (Recommended)
%RHSvector=[];
