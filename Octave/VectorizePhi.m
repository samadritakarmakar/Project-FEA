function phi=VectorizePhi(phiTemp, vectorLevel)
    NoOfPhiColumns=size(phiTemp,2);
    phi=spalloc(vectorLevel,NoOfPhiColumns*vectorLevel,NoOfPhiColumns*vectorLevel);
    for i=1:NoOfPhiColumns
        startColumn=i*vectorLevel-(vectorLevel-1);
        endColumn=startColumn+vectorLevel-1;
        %phi(1:vectorLevel,startColumn:endColumn)=sparse(diag(repmat(phiTemp(i),1,vectorLevel)));
        phi(1:vectorLevel,startColumn:endColumn)=sparse(bsxfun(@times, speye(vectorLevel), phiTemp(i)));
    end
end
