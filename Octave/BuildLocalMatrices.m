function [LHSmatrixLocal RHSmatrixLocal RHSvectorLocal] = BuildLocalMatrices(LocalGaussMatrices, LocalMatrices, F, w)
        LHSmatrixLocal=LocalMatrices.LHSmatrixLocal;
        RHSmatrixLocal=LocalMatrices.RHSmatrixLocal;
        RHSvectorLocal=LocalMatrices.RHSvectorLocal;
        
        LHSmatrixLocalGauss=LocalGaussMatrices.LHSmatrixLocalGauss;
        RHSmatrixLocalGauss=LocalGaussMatrices.RHSmatrixLocalGauss;
        RHSvectorLocalGauss=LocalGaussMatrices.RHSvectorLocalGauss;
        
        
        %Jonathan Whiteley Finite Element Methods A Practical Guide eqn 7.36
        LHSmatrixLocal=w(1)*w(2)*w(3)*LHSmatrixLocalGauss*det(F) + LHSmatrixLocal;
        RHSmatrixLocal=w(1)*w(2)*w(3)*RHSmatrixLocalGauss*det(F) + RHSmatrixLocal;
        RHSvectorLocal=w(1)*w(2)*w(3)*RHSvectorLocalGauss*det(F) + RHSvectorLocal;
        
end
