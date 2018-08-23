load input/node_t;
x=zeros(3*node_t,1);
F=zeros(3*node_t,1);
node_anc=input("How many nodes are anchored?: ");
for i=1:node_anc
 printf("Enter node no. anchored:\n");
 printf("%u: ",i);
 j=input("");
 x(3*j-2)=0;
 x(3*j-1)=0;
 x(3*j)=0;
endfor

node_disp=input("Enter number of displacements known. ");
for i=1:node_disp
 printf("Displacement %u:\n",i);
 chc=input("is it in x direction?(y/n): \n","s");
 if(chc=='y'||chc=='Y')
  j=input("Enter the node no.: ");
  x(3*j-2)=input("Enter displacment for x-axis: ");
 else
  chc=input("is it in y direction?(y/n): \n","s");
  if(chc=='y'||chc=='Y')
   j=input("Enter the node no.: ");
   x(3*j-1)=input("Enter displacment for y-axis: ");
  else
   chc=input("is it in z direction?(y/n): \n","s");
   if(chc=='y'||chc=='Y')
    j=input("Enter the node no.: ");
    x(3*j)=input("Enter displacment for z-axis: ");
   endif
  endif
 endif
endfor

node_frc=input("Enter number of forces known. ");
for i=1:node_frc
printf("Displacement %u:\n",i);
 chc=input("is it in x direction?(y/n): \n","s");
 if(chc=='y'||chc=='Y')
  j=input("Enter the node no.: ");
  F_locatn(i)=3*j-2;
  F(F_locatn(i))=input("Enter force in direction of x-axis: ");
 else
  chc=input("is it in y direction?(y/n): \n","s");
  if(chc=='y'||chc=='Y')
   j=input("Enter the node no.: ");
   F_locatn(i)=3*j-1;
   F(F_locatn(i))=input("Enter force in direction of y-axis: ");
  else
   chc=input("is it in z direction?(y/n): \n","s");
   if(chc=='y'||chc=='Y')
    j=input("Enter the node no.: ");
    F_locatn(i)=3*j;
    F(F_locatn(i))=input("Enter force in direction of z-axis: ");
   endif
  endif
 endif
endfor

save -ascii boundary/F F
save -ascii boundary/x x
save -ascii boundary/F_locatn F_locatn