function ContainsNodes=SortAndKeep(ElementNodes)
    [rows columns]=size(ElementNodes);
    for row=1:rows
        for column=1:columns
            if (row==1 && column==1)
                ContainsNodes(1)=ElementNodes(row,column);
                ContainsNodes(2)=ElementNodes(row,column+1);
            end
            flag=1;
            first=1;
            last=length(ContainsNodes);
            while(flag)
                middle=int32((first+last)/2);
                if (ContainsNodes(middle-1)<ElementNodes(row,column) && ContainsNodes(middle)<ElementNodes(row,column))
                first=middle;
                if (middle==last)
                    ContainsNodes=[ContainsNodes,ElementNodes(row,column)];
                    flag=0;
                end
                elseif (ContainsNodes(middle-1)>ElementNodes(row,column) && ContainsNodes(middle)>ElementNodes(row,column))
                last=middle;
                if (middle-1==first)
                        ContainsNodes=[ElementNodes(row,column),ContainsNodes];
                        flag=0;
                    end
                elseif (ContainsNodes(middle-1)<ElementNodes(row,column) && ContainsNodes(middle)>ElementNodes(row,column))
                ContainsNodes=[ContainsNodes(1:middle-1),ElementNodes(row,column),ContainsNodes(middle:end)];
                flag=0;
                elseif(ContainsNodes(middle-1)==ElementNodes(row,column) || ContainsNodes(middle)==ElementNodes(row,column))
                flag=0;
                end
            end
        end
    end
end 
