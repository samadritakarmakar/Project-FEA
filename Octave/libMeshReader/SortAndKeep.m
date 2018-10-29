function ContainsNodes=SortAndKeep(ElementNodes)
    [rows columns]=size(ElementNodes);
    ContainsNodes=zeros(1,rows+columns);
    %{
    ContainsNodes(1)=ElementNodes(1,1);
    ContainsNodes(2)=ElementNodes(1,2);
    
    for row=1:rows
        for column=1:columns
          %  if (row==1 && column==1)
          %      ContainsNodes(1)=ElementNodes(row,column);
          %      ContainsNodes(2)=ElementNodes(row,column+1);
          %  end
            flag=1;
            first=1;
            last=length(ContainsNodes);
            while(flag)
                middle=int32((first+last)/2);
                %case1=(ContainsNodes(middle-1)<ElementNodes(row,column) && ContainsNodes(middle)<ElementNodes(row,column));
                %case2=(ContainsNodes(middle-1)>ElementNodes(row,column) && ContainsNodes(middle)>ElementNodes(row,column));
                %case3=(ContainsNodes(middle-1)<ElementNodes(row,column) && ContainsNodes(middle)>ElementNodes(row,column));
                %case4=(ContainsNodes(middle-1)==ElementNodes(row,column) || ContainsNodes(middle)==ElementNodes(row,column));
                
                if (ContainsNodes(middle-1)<ElementNodes(row,column) && ContainsNodes(middle)<ElementNodes(row,column))
                    first=middle;
                    if (middle==last)
                        %ContainsNodes=[ContainsNodes,ElementNodes(row,column)];
                        ContainsNodes(1:middle+1)=ElementNodes(row,column);
                        flag=0;
                    end
                elseif (ContainsNodes(middle-1)>ElementNodes(row,column) && ContainsNodes(middle)>ElementNodes(row,column))
                    last=middle;
                    if (middle-1==first)
                        %ContainsNodes=[ElementNodes(row,column),ContainsNodes];
                        ContainsNodes(1:middle+1)=[ElementNodes(row,column),ContainsNodes]
                        flag=0;
                    end
                elseif (ContainsNodes(middle-1)<ElementNodes(row,column) && ContainsNodes(middle)>ElementNodes(row,column))
                    %ContainsNodes=[ContainsNodes(1:middle-1),ElementNodes(row,column),ContainsNodes(middle:end)];
                    ContainsNodes(1:middle+1)=[ContainsNodes(1:middle-1),ElementNodes(row,column),ContainsNodes(middle:end)];
                    flag=0;
                elseif(ContainsNodes(middle-1)==ElementNodes(row,column) || ContainsNodes(middle)==ElementNodes(row,column))
                    flag=0;
                end
            end
        end
    end
%}
for row=1:rows
    ContainsNodes(1,columns*row-(columns-1):columns*row)=ElementNodes(row,:);
end
ContainsNodes=unique(ContainsNodes);
ContainsNodes=ContainsNodes(ContainsNodes~=0);






end 
