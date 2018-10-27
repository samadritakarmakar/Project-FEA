function [LHSmatrix RHSmatrix RHSvector]=UserFunction(x,u,v,gradu,gradv, ElementNum);

LHSmatrix = ScalarDot(gradv,gradu);
RHSmatrix=[];
RHSvector=v.*f(x);
