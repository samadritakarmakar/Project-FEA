function [LHSmatrix RHSmatrix RHSvector]=UserFunction(x,u,v,gradu,gradv, ElementNum);

LHSmatrix = ScalarDot(gradv,gradv);
RHSmatrix=[];
RHSvector=v.*f(x);
