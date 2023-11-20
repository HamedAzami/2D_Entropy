function Out_entropy=EnsDispEn_2D(I,m,nc,MA)
%
% This function calculates ensemble two-dimesnional 
% dispersion entropy (EnsDispEn_2D) and two-dimensional
% dispersion entropy (DispEn_2D) of an image I based on 
% several possible maps
%  
% Inputs:
%
% I: image
% m: embedding dimension
% nc: number of classes (it is usually equal a number to 3 - 9)
% MA: 1: LM   2: NCDF    3: LOGSIG  4: TANSIG   5: SORT
%
% Outputs:
%
% Out_entropy: scalar quantity - the DispEn_2D or EnsDispEn2D of I if MA corresponds to a
% single map or several maps, respectively

% npdf: a vector of length (nc^m)^m, showing the normalized number of disersion patterns of I
%
% Ref:
% [1] H. Azami, L. E. da Silva, A. C. M. Omoto, & A. Humeau-Heurtier, "Two-dimensional dispersion entropy: An information-theoretic
% method  for irregularity analysis of images", Signal Processing: Image Communication, vol. 75, pp. 178-187. 2019.
% [2] , tansig, logsig, and sorting
% [3] Gaudêncio, A. S., Azami, H., Cardoso, J. M., Vaz, P. G., & Humeau-Heurtier, A. (2023). Bidimensional ensemble entropy: Concepts 
% and application to emphysema lung computerized tomography scans. Computer Methods and Programs in Biomedicine, 107855.
%
% If you use the code, please make sure that you cite references [1], [2], and [3].
%
%  19-september-2022
%%
%



Out_LM=NaN;
Out_NCDF=NaN;
Out_LOGSIG=NaN;
Out_TANSIG=NaN;
Out_SORT=NaN;

npdf_LM=NaN*ones(1,(nc^m)^m);
npdf_NCDF=NaN*ones(1,(nc^m)^m);
npdf_LOGSIG=NaN*ones(1,(nc^m)^m);
npdf_TANSIG=NaN*ones(1,(nc^m)^m);
npdf_SORT=NaN*ones(1,(nc^m)^m);


if ismember(1,MA)
    [Out_LM, npdf_LM]=DispEn_2D(I,m,nc,'LM')   ; 
end


if ismember(2,MA)
    [Out_NCDF, npdf_NCDF]=DispEn_2D(I,m,nc,'NCDF')    ;
end


if ismember(3,MA)
    [Out_TANSIG, npdf_TANSIG]=DispEn_2D(I,m,nc,'LOGSIG')   ; 
end


if ismember(4,MA)
    [Out_LOGSIG, npdf_LOGSIG]=DispEn_2D(I,m,nc,'TANSIG')  ;  
end


if ismember(5,MA)
    [Out_SORT, npdf_SORT]=DispEn_2D(I,m,nc,'SORT')    ;
end

npdf=nanmean([npdf_LM;npdf_NCDF;npdf_TANSIG;npdf_LOGSIG;npdf_SORT]);


p=npdf(npdf~=0);
Out_entropy = -sum(p .* log(p));