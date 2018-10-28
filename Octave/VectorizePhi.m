function phi=VectorizePhi(phiTemp, NoOfInterpolatedCoords)
    NoOfPhiColumns=size(phiTemp,2);
    phi=spalloc(NoOfInterpolatedCoords,NoOfPhiColumns*NoOfInterpolatedCoords,NoOfPhiColumns*NoOfInterpolatedCoords);
    for i=1:NoOfPhiColumns
        startColumn=i*NoOfInterpolatedCoords-(NoOfInterpolatedCoords-1);
        endColumn=startColumn+NoOfInterpolatedCoords-1;
        %phi(1:NoOfInterpolatedCoords,startColumn:endColumn)=sparse(diag(repmat(phiTemp(i),1,NoOfInterpolatedCoords)));
        phi(1:NoOfInterpolatedCoords,startColumn:endColumn)=sparse(bsxfun(@times, speye(NoOfInterpolatedCoords), phiTemp(i)));
    end
end
