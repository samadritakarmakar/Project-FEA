function [LHSmatrix RHSmatrix RHSvector]=UserFunction(xInterpolated,VectorizedX,u,v,gradu,gradv, ElementNum);

LHSmatrix = TestFunctionDot(gradv,gradu);
RHSmatrix=[];
RHSvector=VectorDot(v,f(VectorizedX)); %u added to enable integration at Gauss Point (Recommended)
%RHSvector=[];
