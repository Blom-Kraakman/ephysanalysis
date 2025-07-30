function [Vs,Ph,Ray] = calcVS(CycT)
%calcVS Calculate Vector Strength, Phase and Rayleigh Criterion
%   INPUT:
%   CycT = (cell) array of spike time re cycle [0,1)
%
%   OUTPUT:
%   Vs  = scalar; Vector strength
%   Ph  = scalar; Phase in cycles
%   Ray = [in progress] p value of rayleigh criterion
    
    if iscell(CycT) % convert cell array to matrix
        CycT = cell2mat(CycT);
    end
    
    N       =   numel(CycT);
    x		=	mean(cos(2*pi*CycT(:)));
    y		=	mean(sin(2*pi*CycT(:)));
    
    Vs = sqrt(x^2+y^2);
    Ph = mod(atan2(y,x)/2/pi,1);
    Ray = rayleigh_pval(Vs, N);
        
end 

