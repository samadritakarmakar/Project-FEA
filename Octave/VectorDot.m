function product= VectorDot(v,p2)
    %NumOfPhi=size(p2,2);
    %for i=1:NumOfPhi
    %    p1Dash=repmat(p1(:,i),1,NumOfPhi);
    %    productDash=p1Dash.*p2;
    %    product(i,:)=sum(productDash,1);
    %end
    product=v*(v'*p2);

    
end
