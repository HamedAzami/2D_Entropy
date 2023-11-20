function Out = EnsFuzEnM_2D(I, M, r, n)
%
% This function calculates Ensemble 
% Fuzzy Entropy based on multiple embedding dimensions
% (EnsFuzEnM_2D) of an image I 
%  
% Inputs:
% I: Image
% M: embedding dimension vector 
% r: tolerance value
% n:fuzzy power
%
% Outputs:
%
% Out: scalar quantity - the EnsFuzEnM_2D value

% Ref:
% [1] Gaudêncio, A. S., Azami, H., Cardoso, J. M., Vaz, P. G., & Humeau-Heurtier, A. (2023). Bidimensional ensemble entropy: Concepts 
% and application to emphysema lung computerized tomography scans. Computer Methods and Programs in Biomedicine, 107855.
%
% If you use the code, please make sure that you cite references [1].
%
%
% andreia.gaudencio@fis.uc.pt
%  19-september-2022

I = double(I);

Phi1 = zeros(1,length(M));
Phi2 = zeros(1,length(M));

for i_m=1:length(M)
    m = M(i_m);
    
    r = r*std2(I);
    
    [~, Phi] = FuzEn2D(m,r,n,I);
    
    Phi1(i_m) = Phi(1,1);
    Phi2(i_m) = Phi(2,1);
end

Phi1m=mean(Phi1);
Phi2m=mean(Phi2);

Out = log(Phi1m/Phi2m);

end