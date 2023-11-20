function [fuzzy, Phi] = FuzEn2D(m,r,n,I)
%
% This function calculates
% Fuzzy Entropy of an image I 
%  
% Inputs:
% m: embedding dimension 
% r: tolerance value
% n:fuzzy power
% I: Image

%
% Outputs:
%
% fuzzy: scalar quantity of FuzEn2D
% Phi: vector with the number of similarity degree sum for the embedding
% dimension m and m+1

% Ref:
% [1] Gaudêncio, A. S., Azami, H., Cardoso, J. M., Vaz, P. G., & Humeau-Heurtier, A. (2023). Bidimensional ensemble entropy: Concepts 
% and application to emphysema lung computerized tomography scans. Computer Methods and Programs in Biomedicine, 107855.
%
% If you use the code, please make sure that you cite references [1].
%
%
% andreia.gaudencio@fis.uc.pt
%  19-september-2022
%%
[nx, ny]=size(I); %[number of rows, number of columns, number of slices]
Phi=zeros(2,1);

%Total number of patterns (both m and m+1 points)
den=(nx-m)*(ny-m);

templ=zeros(m,m,den);
templ2=zeros(m+1,m+1,den);
cc=0; 

for yi=1:ny-m
for xi=1:nx-m
    
    cc=cc+1;
    
    templ(:,:,cc)=I(xi:xi+m-1,yi:yi+m-1); 
    templ2(:,:,cc)=I(xi:xi+m,yi:yi+m);
    
end 
end 

ss=zeros(1,den);
ss2=zeros(1,den);

parfor jj=1:den-1
    
    
    dist=max(abs(templ(:,:,jj+1:den)-repmat(templ(:,:,jj),[1,1,den-jj])),[],[1 2]);
    dist2=max(abs(templ2(:,:,jj+1:den)-repmat(templ2(:,:,jj),[1,1, den-jj])),[],[1 2]);
    
    D=exp((-dist.^n)./r); 
    D2=exp((-dist2.^n)./r);
    
    ss(jj)=sum(D)/(den-1);
    ss2(jj)=sum(D2)/(den-1);
    
end


Phi(1,1)=sum(ss)/(den);

Phi(2,1)=sum(ss2)/(den);


fuzzy=log(Phi(1,1)/Phi(2,1)); 

end