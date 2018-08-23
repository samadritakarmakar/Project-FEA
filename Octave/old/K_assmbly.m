load input/elem_t;
load input/node_t;
load input/node1;
load input/node2;
load input/type;
load input/theta_x;
load input/theta_y;
load input/theta_z;
load input/A;
load input/E;
load input/L;
load input/k_stiff;

K=zeros(3*node_t);
for i=1:elem_t
 if(type(i)==1)
  K=K+k_spring(k_stiff(i),theta_x(i),theta_y(i),theta_z(i),node1(i),node2(i),node_t);
 elseif (type(i)==2)
  K=K+k_truss(A(i),E(i),L(i),theta_x(i),theta_y(i),theta_z(i),node1(i),node2(i),node_t);
 endif
endfor
save -ascii boundary/K K
