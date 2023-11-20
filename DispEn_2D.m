function [Out_DisEn npdf]=DispEn_2D(I,m,nc,MA)
%
% This function calculates two--dimesnional dispersion entropy (DisEn_2D) of an image I
%
% Inputs:
%
% I: image
% m: embedding dimension
% nc: number of classes (it is usually equal a number to 3 - 9)

%
% Outputs:
%
% Out_DisEn: scalar quantity - the DisEn_2D of I
% npdf: a vector of length (nc^m)^m, showing the normalized number of disersion patterns of I
%
% Ref:
% [1] H. Azami, L. E. da Silva, A. C. M. Omoto, & A. Humeau-Heurtier, "Two-dimensional dispersion entropy: An information-theoretic
% method  for irregularity analysis of images", Signal Processing: Image Communication, vol. 75, pp. 178-187. 2019.
% [2] M. Rostaghi and H. Azami, "Dispersion Entropy: A Measure for Time-Series Analysis", IEEE Signal Processing Letters. vol. 23, n. 5, pp. 610-614, 2016.
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Anne Humeau-Heurtier and Hamed Azami
% Emails: anne.humeau@univ-angers.fr and hmd.azami@gmail.com
%
%  27-march-2019
%%
%


I=double(I);
N_x=size(I,1);
N_y=size(I,2);

sigma_I=std2(I);
mu_I=mean2(I);


switch MA
    case   'LM'
        y=mapminmax(I(:)',0,1);
        y(y==1)=1-eps;
        y(y==0)=eps;
        y=reshape(y,N_x,N_y);
        z=round(y*nc+0.5);
        
        
    case 'NCDF'
        y=normcdf(I(:)',mu_I,sigma_I);
        y=mapminmax(y,0,1);
        y(y==1)=1-eps;
        y(y==0)=eps;
        y=reshape(y,N_x,N_y);
        z=round(y*nc+0.5);
        
    case 'LOGSIG'
         y=logsig((I(:)'-mu_I)/sigma_I);
        y=mapminmax(y,0,1);
        y(y==1)=1-eps;
        y(y==0)=eps;
        y=reshape(y,N_x,N_y);
        z=round(y*nc+0.5);
        
    case 'TANSIG'
        
         y=tansig((I(:)'-mu_I)/sigma_I);
        y=mapminmax(y,0,1);
        y(y==1)=1-eps;
        y(y==0)=eps;
        y=reshape(y,N_x,N_y);
        z=round(y*nc+0.5);
        
    case 'SORT'
        if size(I,1)>size(I,2)
           I=I'; 
        end
        
        for ii=1:size(I,1)
            x=squeeze(I(ii,:));
       N=length(x);
       x=x(1:nc*floor(N/nc));
       N=length(x); 
       [sx osx]=sort(x);         
      
        Fl_NC=N/nc;
        cx=[];
        for i=1:nc
            cx=[cx i*ones(1,Fl_NC)];
        end
        for i=1:N
            z(ii,i)=cx(osx==i);
        end
        end
        if size(I,1)>size(I,2)
           z=z'; 
        end
        
 end



N_x=size(z,1);
N_y=size(z,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if m = 1
if(m==1)
    % generation of the matrices containing all possible submatrices
    all_patterns=NaN*ones(1,(nc^m)^m);
    index_all_patterns=1;
    for i = 1:nc^m,
        all_patterns(index_all_patterns) = i;
        index_all_patterns = index_all_patterns+1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if m = 2:
elseif(m==2)
    
    subpatterns=zeros(1,nc^m);
    index=1;
    %generation of all possible sub_patterns (of size m)
    for i = 1:nc, % m loops -> one on i, one on j, and one on k
        for j = 1:nc,
            subpatterns(index)=10*i+j;
            index=index+1;
        end
    end
    
    % generation of the matrices containing all possible submatrices
    all_patterns=NaN*ones(1,(nc^m)^m);
    index_all_patterns=1;
    for i = 1:nc^m,
        for j = 1:nc^m,
            all_patterns(index_all_patterns) = 100*subpatterns(i)+subpatterns(j);
            index_all_patterns = index_all_patterns+1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % You can add any other m-values by adapting the following code :
elseif(m==3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if m = 3:
    
    subpatterns=zeros(1,nc^m);
    index=1;
    %generation of all possible sub_patterns (of size m)
    for i = 1:nc, % m loops -> one on i, one on j, and one on k
        for j = 1:nc,
            for k = 1:nc,
                subpatterns(index)=100*i+10*j+k;
                index=index+1;
            end
        end
    end
    
    % generation of the matrices containing all possible submatrices
    all_patterns=NaN*ones(1,(nc^m)^m);
    index_all_patterns=1;
    for i = 1:nc^m,
        for j = 1:nc^m,
            for k = 1:nc^m,
                all_patterns(index_all_patterns) = 1000000*subpatterns(i)+1000*subpatterns(j)+subpatterns(k);
                index_all_patterns = index_all_patterns+1;
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
elseif(m==4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if m = 4:
    
    subpatterns=zeros(1,nc^m);
    index=1;
    %generation of all possible sub_patterns (of size m)
    for i = 1:nc, % m loops -> one on i, one on j, one on k, and one on h
        for j = 1:nc,
            for k = 1:nc,
                for h = 1:nc,
                    subpatterns(index)=1000*i+100*j+10*k+h;
                    index=index+1;
                end
            end
        end
    end
    
    % generation of the matrices containing all possible submatrices
    all_patterns=ones(1,(nc^m)^m);
    index_all_patterns=1;
    for i = 1:nc^m,
        for j = 1:nc^m,
            for k = 1:nc^m,
                for h = 1:nc^m,
                    all_patterns(index_all_patterns) = 10^12*subpatterns(i)+10^8*subpatterns(j)+10^4*subpatterns(k)+subpatterns(h);
                    index_all_patterns = index_all_patterns+1;
                end
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif(m==5)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if m = 5:
    
    subpatterns=zeros(1,nc^m);
    index=1;
    %generation of all possible sub_patterns (of size m)
    for i = 1:nc,
        for j = 1:nc,
            for k = 1:nc,
                for h = 1:nc,
                    for g = 1:nc,
                        subpatterns(index)=10000*i+1000*j+100*k+10*h+g;
                        index=index+1;
                    end
                end
            end
        end
    end
    
    % generation of the matrices containing all possible submatrices
    all_patterns=NaN*ones(1,(nc^m)^m);
    index_all_patterns=1;
    for index_all_patterns_i = 1:nc^m,
        for index_all_patterns_j = 1:nc^m,
            for index_all_patterns_k = 1:nc^m,
                for index_all_patterns_h = 1:nc^m,
                    for index_all_patterns_g = 1:nc^m,
                        all_patterns(index_all_patterns) = 10^20*subpatterns(i)+10^15*subpatterns(j)+10^10*subpatterns(k)+10^5*subpatterns(h)+subpatterns(g);
                        index_all_patterns = index_all_patterns+1;
                    end
                end
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


% for all m values:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npdf=zeros(1,(nc^m)^m); %initialize frequency array

ind=1;

% We extract all dispesion patterns of I
for i_z=1:(N_x-m+1)
    for j_z=1:(N_y-m+1)
        extract_z=z(i_z:i_z+m-1,j_z:j_z+m-1);
        Ext_z(ind)=extract_z(1);
        % We map each m*m matrix to a number. For rxample [1 2; 3 4] is mapped to 1324.
        for i_dm=1:m^2-1
            Ext_z(ind)=extract_z(i_dm+1)*10^i_dm+Ext_z(ind);
        end
        ind=ind+1;
    end
end


for i_f=1:length(all_patterns)
    [R,C]=find(Ext_z==all_patterns(i_f));
    npdf(i_f)=length(R)/((N_x-m+1)*(N_y-m+1));
end

p=npdf(npdf~=0);
Out_DisEn = -sum(p .* log(p));
