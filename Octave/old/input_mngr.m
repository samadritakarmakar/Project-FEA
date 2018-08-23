i=1;
node_t=0;
elem_t=0;
do
 type(i)=input("Enter\n 1: for Spring\n 2: for Truss: ");
 if(type(i)==1)
  A(i)=0;
  E(i)=0;
  L(i)=0;
  k_stiff(i)=input("Enter Spring Stiffness: ");
  theta_x(i)=input("Angle with x-axis: ");
  theta_y(i)=input("Angle with y-axis: ");
  theta_z(i)=input("Angle with z-axis: ");

 
 elseif(type(i)==2)
  A(i)=input("Area : ");
  E(i)=input("Elasticity: ");
  k_stiff(i)=0;
  ch=input("Enter\n 1: for coordinate entry\n 2: for angle and length entry: ");
  if(ch==1)
   x1=("x1: ");
   y1=("y1: ");
   z1=("z1: ");
   x2=("x2: ");
   y2=("y2: ");
   z2=("z2: ");
   L(i)=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
   theta_x(i)=acosd((x2-x1)/L);
   theta_y(i)=acosd((y2-y1)/L);
   theta_z(i)=acosd((z2-z1)/L);
  else
   L(i)=input("Enter Length: ");
   theta_x(i)=input("Angle with x-axis: ");
   theta_y(i)=input("Angle with y-axis: ");
   theta_z(i)=input("Angle with z-axis: ");
  endif
 endif
 
 node1(i)=input("Node 1: ");
 node2(i)=input("Node 2: "); 
  if(node1(i)>node_t)
  node_t=node1(i);
 elseif (node2(i)>node_t)
  node_t=node2(i);
 endif
 elem_t=i;
 i++;
 ch=input("Enter\n 0: to Stop entry\n 1: to enter next element: ");
until (ch==0)

save -ascii input/elem_t elem_t;
save -ascii input/node_t node_t;
save -ascii input/node1 node1;
save -ascii input/node2 node2;
save -ascii input/type type;
save -ascii  input/theta_x theta_x;
save -ascii input/theta_y theta_y;
save -ascii input/theta_z theta_z;
save -ascii input/A A;
save -ascii input/E E;
save -ascii input/L L;
save -ascii input/k_stiff k_stiff;
 