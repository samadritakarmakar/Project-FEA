%Used to put values of stiffness of Truss in the global Stiffness matrix K.
%Takes Area, Elasticity, Length, angle with x axis, node number 1, node number 2 and 
%total number of nodes as argument.
function k=k_truss(A,E,L,theta_x,theta_y,theta_z,node1,node2,node_t)
Cx=cosd(theta_x);
Cy=cosd(theta_y);
Cz=cosd(theta_z);
k=zeros(3*node_t);
%The right positions of k generally ranges from 3n-2 to 3n where n is node number.
n1=3*node1;
n2=3*node2;
lambda=A*E/L*[Cx^2, Cx*Cy, Cx*Cz;
              Cy*Cx, Cy^2, Cy*Cz;
              Cz*Cx, Cz*Cy, Cz^2];
%Puts the values in the right positions.
k(n1-2:n1,n1-2:n1)=lambda;
k(n1-2:n1,n2-2:n2)=-lambda;
k(n2-2:n2,n1-2:n1)=-lambda;
k(n2-2:n2,n2-2:n2)=lambda;
endfunction
