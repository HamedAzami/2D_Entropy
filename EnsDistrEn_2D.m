function Out= EnsDistrEn_2D(I,M,B)

%
% This function calculates Ensemble 
% Distribution Entropy (EnsDistEn_2D)
% of an image I 
%  
% Inputs:
% I: Image
% M: embedding dimension vector 
% B: Histogram has B bins. B is an intiger power of 2 
%
% Outputs:
%
% Out: scalar quantity - the EnsDistEn_2D value

% Ref:
% [1] Gaudêncio, A. S., Azami, H., Cardoso, J. M., Vaz, P. G., & Humeau-Heurtier, A. (2023). Bidimensional ensemble entropy: Concepts 
% and application to emphysema lung computerized tomography scans. Computer Methods and Programs in Biomedicine, 107855.
%
% If you use the code, please make sure that you cite references [1].
%
%  19-september-2022

I=double(I); %alteração minha
for i=1:length(M)
    m=M(i);
    Out(i) = DistEn_2D(I,m,B);
end
Out=nanmean(Out);

