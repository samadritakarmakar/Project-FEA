 function j=jacobian(f,x1,OtherData) %pointer of function and column of vectors
 x=vec(x1); %Ensures a column vector
 %x_h=.000000002; %column of vectors  x_h
 x_h_alt=.000000002; %in case and element of x is zero
 x_h=sqrt(eps(x)).*x; %find x_h
 positionOf0=find(~x_h); %find position of x=0
 if ~isempty(positionOf0) %check if none are zero 
    x_h(positionOf0)=eps(1); %set 0 to some value
 end
 l=length(x); %length of x must be equal to length of f(x)
 xPlus_h=repmat(x,1,l)+diag(x_h);

    for i=1:l
        fx1=(feval(f,xPlus_h(:,i),OtherData)-feval(f,x,OtherData));
        fx=vec(fx1); %Ensures a column vector
        j(:,i)=fx/x_h(i); %(f(x_i+ xh_i) - f(x_i))/x_h(i)
    end
 end
