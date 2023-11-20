function [entropy,A,B,N] = EnsSampEn_2D(X,m,R,varargin)
%ensemble SAMPEN2D Two-dimensional Sample Entropy (SampEn2D).
%
%    E = SAMPEN2D(X,M,R) calculates the sample entropy 2D for a
%    bidimensional matrix X, as representing a grayscale image. X will be
%    converted to double type. Calculation is performed using parameters M,
%    the length of the pattern (i.e., squared windows of size M), and R,
%    the matching tolerance factor (see Name-Value 'RPercentage' below).
%
%
%    E = SAMPEN2D(...,D) uses a spatial delay D between pattern elements.
%    Default is 1, so SAMPEN2D(X,M,R) is the same as SAMPEN2D(X,M,R,1).
%
%
%    E = SAMPEN2D(...,'RPercentage',BOOL) controls how the value of R is
%    considered in tolerance factor determination. BOOL is a logical value:
%
%        true    - (default) R value is considered as a percentage. Then,
%                  the tolerance is calculated as R*std2(X).
%        false   - R is considered an absolute value and is directly used
%                  as the tolerance.
%
%    E = SAMPEN2D(...,'MexFunction',BOOL) defines if matching calculation
%    is done by executing a mex function or not. BOOL is a logical value:
%        
%        true    - (default) mex function is executed (faster option). If 
%                  the mex file is not found in the path, local matlab 
%                  function is executed instead.
%        false   - local matlab function is executed (slower option).                     
%
%
%    [E,A,B,N] = SAMPEN2D(X,M,R,...) returns the sample entropy 2D and also
%    the number of counts for matches of patterns of size M (B) and matches 
%    of patterns of size M+1 (A), and the total number of paired patterns N.
%
%
%  
%  This implementation is based on the algorithm published in [1].
%
%  If you use this code preparing a paper for publication, please add a 
%  reference to the following:
%
%  [1] LEV Silva, ACS Senra Filho, et al., "Two-dimensional sample entropy:
%      assessing image texture through irregularity", Biomedical Physics & 
%      Engineering Express, 2(4), p.045002, 2016. 
%      DOI: 10.1088/2057-1976/2/4/045002   
%

% Copyright 2016-2018 Juliano J. Duque
% Licensed under the MIT License
%

    % Parsing input arguments
    % Parsing input arguments
    p=inputParser;

    defaultDelay = 1;
    defaultRPercentage = true;
    defaultMexFunction = true;

    addRequired(p,'X',@isnumeric);
    addRequired(p,'m',@(x) rem(x,1)==0 && (x>0)); %positive integer
    addOptional(p,'d',defaultDelay, @(x) rem(x,1)==0 && (x>0));

    addParameter(p,'RPercentage',defaultRPercentage,@islogical);
    addParameter(p,'MexFunction',defaultMexFunction,@islogical);

    
    for i_r=1:length(R)
    r=R(i_r);
    if i_r==1
    addRequired(p,'r',@(x) isreal(x) && (x>0)); %positive real
    end
    parse(p,X,m,r,varargin{:});
    % Parsed delay value
    d = p.Results.d;
    
    % Convert image to double
    X = double(X);
    
    % Check tolerance factor
    if p.Results.RPercentage
      r = r * std2(X); % Update r
    end

    % Check if mex file exists
    runMex = p.Results.MexFunction;
    if runMex
        if exist('cmatches2d','file')~=3 % mex-file not found in path
            warning('backtrace','off');
            warning('MEX-file "cmatches2d" not found in MATLAB search path.');
            warning('backtrace','on');
            disp('Calculating using local function "matches2d" instead...');
            runMex = false;
        end
    end
        
    % Run pattern matching
    if runMex
        [A(i_r),B(i_r),N] = cmatches2d(X,m,r,d);
    else
        [A(i_r),B(i_r),N] = matches2d(X,m,r,d);
    end

    % Entropy value

    end
    entropy = -log(double(A)/double(B));


end

% Local function that performs the matching
function [A,B,N] = matches2d(I,m,r,d)
    ni = size(I,1); % Number of rows
    nj = size(I,2); % Number of columns

    % Matches counts
    B = uint64(0); % for m
    A = uint64(0); % for m+1

    % Check m and d combination
    if m*d >= min(ni,nj) || m*d+d >= max(ni,nj) 
        %Shortest dimension cannot be less than or equal to m*d (no pattern building)
        %Longest dimension cannot be less than or equal to m*d+d, (no comparisons can be done)
        B = NaN;
        A = NaN;
        N = 0;
        return; 
    end

    % Limit indexes for cols and rows
    ilim = ni - m*d;
    jlim = nj - m*d;
    
    % Total number of patterns (for both m and m+1)
    Nm = ilim * jlim;

    % Total number of comparisons (paired patterns)
    Mh=min(d-1,ilim);
    Mw=min(d-1,jlim);
    N= ((Nm*(Nm-1))/2) - ( ...
         (Nm*2*d*(d-1)) - ...
         ((2*d-1)/2)*( (ilim*Mw*(2*d-1-Mw))+(jlim*Mh*(2*d-1-Mh)) ) + ...
         Mh*Mw*(2*d-1-Mh)*(2*d-1-Mw)/2 ...
       );
    %Note: the equation below holds for square images (ni==nj)
    %N2 = ((Nm*(Nm-1))/2)-(d*(d-1)*((2*Nm)-((2*d-1)/2)*(ilim+jlim)+(d*(d-1)/2)));

    % Pre-allocating index shifts for patterns matching
    n_m = m*m;    % number of pixels m length pattern
    n_m1 = 2*m+1; % number of extra pixels to complete m+1 length pattern
    
    % arrays with shifts to be applied for m and m+1 matching
    s_m = zeros(n_m,1);
    s_m1 = zeros(n_m1,1);
    
    k=1;
    for i = 1:m
        kk = ni * (i-1) * d;
        for j = 1:m
            s_m(k) = kk + (j-1)*d;
            k=k+1;
        end
    end
    k=1;
    for j = 1:m+1 % Compares column m+1      
        s_m1(k) = (ni*m +(j-1)) * d;
        k=k+1;
    end
    for i = 1:m  % Compares row m+1
        s_m1(k) = (ni*(i-1) + m) * d;
        k=k+1;
    end
    % End pre-allocating

    % Note: Matlab sees a 2D array (matrix) as a long single column array.

    % Starts running
    for j1 = 0:jlim-1 % for each column
        p1lim = (j1*ni)+ni-(m*d); % set limit index for template pattern            
        
        p1 = j1*ni+1; % p1 is the template pattern
        i1 = 1; %row iterator
        
        while p1 <= p1lim
            % p2 is the moving pattern
            p2 = p1+d; % constraint to avoid spatial correlation effect      
            i2 = i1+d; % row iterator
            
            for j2 = j1:jlim-1 % for each column
                p2lim = (j2*ni)+ni-(m*d); % set limit index for moving pattern
             
                while p2 <= p2lim
                    % Constraint to avoid spatial correlation effect
                    if (j2-j1) < d && abs(i2-i1) < d  
                        i2=i2+1;
                        p2=p2+1;                    
                        continue;
                    end
                    
                    %pattern matching for m with tolerance r
                    for k=1:n_m
                        dif = abs(I(p1+s_m(k)) - I(p2+s_m(k)));
                        if dif >= r
                            dif = -1;
                            break;
                        end
                    end

                    if dif ~= -1 % are patterns similar for m?
                        B=B+1;
                        
                        % pattern matching for m+1 with tolerance r
                        for k=1:n_m1
                            dif = abs(I(p1+s_m1(k)) - I(p2+s_m1(k)));
                            if dif >= r
                                dif = -1;
                                break;
                            end
                        end
                        
                        if dif ~= -1 % are they similar for m+1?
                            A=A+1;
                        end
                    end

                    i2=i2+1;
                    p2=p2+1;
                end % while p2

                i2=1;
                p2=(j2+1)*ni+1;
            end % for j2
            
            i1=i1+1;
            p1=p1+1;
        end % while p1
    end % for j1

    % Average probabilities Um and Um+1 may be obtained
    % dividing B and A by N, respectively.

end




    