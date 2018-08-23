function fx=f(xn,t)
%fx=xn.^2;

x=xn(1);
y=xn(2);
z=xn(3);

fx=[x^2 + x*y*z + y^2*z^3 - 3; ...
    x*y^2 - y*z^2 - 2*x^2; ...
    x^2*y + y^2*z + z^4 - 4];

    end
