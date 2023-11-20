function Out = EnsPerEn_2D(I,m_min,m_max,tau)
% 
% Calculate the Ensemble 2D permutation entropy (EnsPerEn2D)
% Input:    I: image;
%           m_min: minimum order of 2D permutation entropy
%           m_max: maximum order of 2D permutation entropy
%           tau: delay time of permutation entropy,
% Output:
%           Out: EnsPerEn2D
%

% Ref:
% [1] Gaudêncio, A. S., Azami, H., Cardoso, J. M., Vaz, P. G., & Humeau-Heurtier, A. (2023). Bidimensional ensemble entropy: Concepts 
% and application to emphysema lung computerized tomography scans. Computer Methods and Programs in Biomedicine, 107855.
%
% If you use the code, please make sure that you cite references [1].
%
% andreia.gaudencio@fis.uc.pt
%  19-september-2022

I=double(I);
i_ind=1;
for i=m_max:-1:m_min  
    Out(i_ind) = Perm_entropy_2D(I,i,i,tau,tau);

i_ind=i_ind+1;
end
Out=nanmean(Out);

