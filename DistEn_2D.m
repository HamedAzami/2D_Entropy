function Out_DistEn_2D = DistEn_2D(image,m,B)
%
% This function calculates
% Distribution Entropy (DistEn_2D)
% of an image I 
%  
% Inputs:
% I: Image
% M: embedding dimension vector 
% B: Histogram has B bins. B is an intiger power of 2 
%
% Outputs:
%
% Out: scalar quantity - the DistEn_2D value

% Ref:
% [1] Gaudêncio, A. S., Azami, H., Cardoso, J. M., Vaz, P. G., & Humeau-Heurtier, A. (2023). Bidimensional ensemble entropy: Concepts 
% and application to emphysema lung computerized tomography scans. Computer Methods and Programs in Biomedicine, 107855.
%
% If you use the code, please make sure that you cite references [1].
% 
% andreia.gaudencio@fis.uc.pt
% hmd.azami@gmail.com
%  19-september-2022

image = (image - min(image(:))) ./ (max(image(:)) - min(image(:)));

N_x = size(image,1); % Number of rows
N_y = size(image,2); % Number of columns

for i_z=1:(N_x-m+1)
    for j_z=1:(N_y-m+1)
        temp=image(i_z:i_z+m-1,j_z:j_z+m-1);
        extract_z(i_z,j_z,:)=temp(:);
    end
end

extract_z_2=reshape(extract_z,(N_x-m+1)*(N_y-m+1),m^2);
dv  = pdist(extract_z_2, 'chebychev');
num  = hist(dv, linspace(0, 1, B));


freq = num./sum(num);
p=freq(freq~=0);
Out_DistEn_2D= -sum(p.*log2(p)) ./ log2(B);




