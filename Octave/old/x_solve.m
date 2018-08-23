%Solves the values of column matrix 'x' in the Equation [F]=[k]*[x]
%Takes present x matrix, positions of known F values in [F] matrix,
%global stiffness matrix 'K' and Force matrix [F] as arguments
function x=x_solve(x,F_locatn,K,F)
 %finds the total number of elements in F_locatn
 l=numel(F_locatn);
 %[F]=[F]-[K]*[x]. Equivalent to substitution of known x values and
 %subtraction from Range.
 F=F-K*x;
 K_stub=zeros(l);
 F_stub=zeros(l,1);
 x_stub=zeros(l,1);
 for i=1:l
 %Copies known [F] values to a smaller matrix
 F_stub(i,1)=F(F_locatn(i));
  for j=1:l
   %Copies [K] values with row and column number same as
   %known [F] row numbers to a smaller matrix.
   %used to find unknown values of [x]
   K_stub(i,j)=K(F_locatn(i),F_locatn(j));
  end
 end
 %Equivalent to 'unknown [x]'='[K]^(-1) * '[F] known' 
 x_stub=K_stub\F_stub;
 
 for i=1:l
  %Evaluated 'unknown [x]' values put back to place.
  x(F_locatn(i))=x_stub(i);
 end

   
