function NodePositons=GetNodePostions(ElementNodes,vectorLevel)
    row=size(ElementNodes,1);
    column=size(ElementNodes,2);
    NodePositons=zeros(row,vectorLevel*column);
    for i=1:column
        for j=0:vectorLevel-1
            NodePositons(:,vectorLevel*i-j)=vectorLevel*ElementNodes(:,i)-j*ones(size(row,2),1);
        end
    end
end
