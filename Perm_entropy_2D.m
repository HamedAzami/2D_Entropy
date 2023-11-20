function Out=Perm_entropy_2D(I,wl_x,wl_y,lag_x,lag_y)
% 
% Calculate the Ensemble 2D permutation entropy (EnsPerEn2D)
% Input:    I: image;
%           wl_x: window defined in the x-axis to 
%                 obtain the final embedding dimension (m_x)
%           wl_y: window defined in the y-axis to 
%                 obtain the final embedding dimension (m_y)
%      
%           lag_x: delay time of permutation entropy in the x-axis
%           lag_y: delay time of permutation entropy in the y-axis
% Output:
%           Out: Perm_entropy_2D
%
% Ref:
% [1] Gaudêncio, A. S., Azami, H., Cardoso, J. M., Vaz, P. G., & Humeau-Heurtier, A. (2023). Bidimensional ensemble entropy: Concepts 
% and application to emphysema lung computerized tomography scans. Computer Methods and Programs in Biomedicine, 107855.
%
% If you use the code, please make sure that you cite references [1].
%
% andreia.gaudencio@fis.uc.pt
% hmd.azami@gmail.com
%  19-september-2022


[l_y,l_x]=size(I);
m_x=l_x-(wl_x-1)*lag_x;
m_y=l_y-(wl_y-1)*lag_y;
iv=zeros(wl_x*wl_y,m_x,m_y);
D_emb=wl_x*wl_y;

for j=1:1:m_y;
     for k=1:1:m_x;
     [~,iv(:,k,j)]=sort(reshape((I(j:lag_y:j+(wl_y-1)*lag_y,k:lag_x:k+(wl_x-1)*lag_x)'),[D_emb,1]));
     end
end

iv=reshape(iv,[D_emb,m_x*m_y]);
[~,~,indcs]=unique(iv','rows');
hist_indcs=hist(indcs,1:1:factorial(D_emb));
pk=hist_indcs/(m_x*m_y);
pk=pk(pk~=0);

Out=-dot(pk,log(pk))/log(factorial(D_emb));