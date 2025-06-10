function [X,Res] = SUnWNS_TV(A,Y,varargin)

%% [X,Res] = SUnWNS_TV(A,Y,varargin)
%  
%  SUnWNS_TV -> Stain Unmixing with nonnegativity, Weight Nucleus Sparsity 
%  and Total Variation is a stain unmixing method proposed for Papanicolaou 
%  stain estimation from RGB images, which is introduced in
% 
%  N. Gong et al., “Papanicolaou Stain Unmixing for RGB Image Using 
%  Weighted Nucleus Sparsity and Total Variation Regularization”.
%  
%  The MATLAB codes were modified based on SUNSAL-TV proposed by 
%  M.-D. Iordache, J. Bioucas-Dias, and A. Plaza (2012).
%
%
%% ----------------------------- Description ------------------------------
%
%  Papanicolaou stain discriminates cellular components by staining nuclei 
%  with hematoxylin (H), cytoplasm with Eosin Y (EY), Light Green SF 
%  yellowish (LG), and Orange G (OG) and lipids with Bismarck brown Y (BY).
%  The number of dyes surpasses the number of RGB channels. We propose an 
%  optimization-based stain unmixing method for Papanicolaou-stained RGB 
%  images, incorporating three priors of stain abundance:
%
%  Nonnegativity: At any given pixel, stain abundance and its corresponding
%                 OD value must be nonnegative.
%
%  Weighted Nucleus Sparsity: Hematoxylin (H) typically binds to nuclei, 
%                             whereas the cytoplasm is nearly unstained by 
%                             H. To enhance the sparsity of H, we introduce
%                             a weight that iteratively assigns higher 
%                             constraints to pixels with lower H abundances
%                             and lower constraints to higher H abundances.    
%
%  Piecewise Smoothness: Stain abundances are expected to vary smoothly 
%                        within homogeneous regions.
%
%  
%  SUnWNS_TV solves the following optimization problem:
%
%     Definitions:
%
%      A  -> 3 * r: (3: number of RGB channels, r: number of dyes) 
%                   stain absorption coefficient matrix; each row
%                   represents the normalized Red/Green/Blue absorbance of
%                   one dye.
%
%      X  -> r * N: (r: number of dyes, N: number of pixels)
%                   stain abundance matrix; each row represents the 
%                   1*N abundance vector of a correspondent dye 
%                   (1-4: EY, H, LG and OG), each column contains the 
%                   abudance of a correspondent pixel.
%
%      Y  -> 3 * N: absorbance matrix 
%                   (optical density computed from RGB image using 
%                   Beer-Lambert Law)
%
%      Note: We excluded BY from the study and considered only the effects of 
%      the remaining four dyes (EY, H, LG, and OG) since BY is almost 
%      transparent, resulting in r = 4.
%
%      Optimization problem:
%
%      min  (1/2) ||A X-Y||^2_F  + i_R_+(X) + lambda_1  w ⊙ ||Xh||_1
%       X                         + lambda_tv ||HX||_{1,1};
%
%
%      where
%
%        Xh -> 1 * N: Xh = X(2,:)
%                     stain abundance vector of Hematoxylin
%
%        w ->  1 * N: w_t+1 = exp(Xh_t)
%                     weights for nucleus spasity. The nucleus sparsity 
%                     weight vector at iteration 𝑡 + 1 is computed from
%                     the stain abundance vector of H at iteration t.
%
%        ⊙: elementwise product
%
%        (1/2) ||A X-Y||^2_F is a quadratic data misfit term
%
%        i_R_+(𝐗) = Σi_R_+(x_i) is the indicator function of the nonnegative 
%        real number set R_+ (i_R_+(x_i) is zero if 𝐱_i belongs to the 
%        nonnegative orthant and +∞ otherwise).
%
%        ||Xh||_1 is the standard l1 regularizer
%
%        ||HX||_{1,1} is the TV (isotropic regularizer)
%
%
%         H is a linear operator that computes the horizontal and the
%         vertical differences on each band of X.  
%         Let Hh: R^{n_lin*n_col}-> R^{n_lin*n_col} be a linear operator that computes the 
%         horizontal first order differences per band. HhX computes a 
%         matrix of the same size of X, where [HhX](i,j) = X(i,h(j))-X(i,j),
%         where h(j) is the index of pixel on the right hand side of j.
%         For the vertical differnces, we have a similar action of Lv:
%         [HvX](i,j) = X(v(i),j)-X(i,j), where v(i) is the index of pixel
%         on the top hand side of j.
%         We consider the anisotropic type of Total variation:
%         ||HX||_{1,1} := ||[Hh; Hv]X||_{1,1}
%
%
%
% -------------------------------------------------------------------------
%
%
%
%% -------------------- Line of Attack  -----------------------------------
%
%  SUnWNS_TV solves the above optimization problem by introducing a variable
%  splitting and then solving the resulting constrained optimization with
%  the augmented Lagrangian method.
%
%
%   The initial problem is converted into
%
%    min  (1/2) ||A X-Y||^2_F  + i_R_+(X)
%     X                        + lambda_1  W ⊙ ||X||_{1,1}
%                              + lambda_tv ||HX||_{1,1};
%
%
%  Then, we apply the following variable splitting
%
%
%    min  (1/2) ||V1-Y||^2     + lambda_1  W ⊙ ||V2||_{1,1}
%  X,V1,V2,V3,V4,V5            + lambda_tv ||V4||_{1,1}
%                              + i_R_+(V5)
%                              
%
%     subject to:  V1  =  AX
%                  V2  =  X
%                  V3  =  X
%                  V4  =  HV3
%                  V5  =  X
%
%
%  For details, please refer to
%
%  N. Gong et al., “Papanicolaou Stain Unmixing for RGB Image Using 
%  Weighted Nucleus Sparsity and Total Variation Regularization”.
%
%
%
%
% ------------------------------------------------------------------------
%%  ===== Required inputs =============
%
%  A - [3(no. of RGB channels) * 4 (no. of dyes)] 
%      stain absorption coefficient matrix.
%
%  Y - [3(no. of RGB channels) x N(pixels)]
%      optical density.
%
%
%%  ====================== Optional inputs =============================
%
%
%  'LAMBDA_1' - regularization parameter for l11 norm.
%               Default: 0;
%
%  'LAMBDA_TV' - regularization parameter for TV norm.
%                Default: 0;
%
%  'IM_SIZE'   - [nlins, ncols]   number of lines and rows of the
%                observation.
%                Note:  n_lin*n_col = N
%
%  'AL_ITERS' - (double):  Maximal number of augmented Lagrangian iterations
%                           Default 2000;
%
%  'MU' - (double):   augmented Lagrangian weight
%                           Default 0.01;
%
%  'POSITIVITY'  = {'yes', 'no'}; Default 'yes'
%                  Enforces the positivity constraint: x >= 0
%
%  'VERBOSE'   = {'yes', 'no'}; Default 'no'
%                 'no' - work silently
%                 'yes' - display warnings
%
%%  =========================== Outputs ==================================
%
%  X  =  [4xn] estimated  X  matrix
%
%

%%
%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
if (nargin-length(varargin)) ~= 2
    error('Wrong number of required parameters');
end
% mixing matrix size
[LM,r] = size(A);
% data set size
[L,N] = size(Y);
if (LM ~= L)
    error('mixing matrix A and data set Y are inconsistent');
end

%%
%--------------------------------------------------------------
% Set the defaults for the optional parameters
%--------------------------------------------------------------
%

% 'LAMBDA_1'
%  l1 regularization
reg_wns = 0;

% 'LAMBDA_TV'
%  TV regularization
reg_TV = 0;
im_size = []; % image size

% 'AL:ITERS'
% maximum number of AL iteration
AL_iters = 2000;

% 'MU'
mu = 0.01;

% 'VERBOSE'
verbose = 'off';

% 'POSITIVITY'
% Positivity constraint
reg_pos = 0;

% initialization
X0 = 0;
Res = [];

%%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------


%--------------------------------------------------------------
% Read the optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'LAMBDA_1'
                lambda_1 = varargin{i+1};
                if lambda_1 < 0
                    error('lambda must be positive');
                elseif lambda_1 > 0
                    reg_wns = 1;
                end
            case 'LAMBDA_TV'
                lambda_TV = varargin{i+1};
                if lambda_TV < 0
                    error('lambda must be non-negative');
                elseif lambda_TV > 0
                    reg_TV = 1;
                end
            case 'IM_SIZE'
                im_size = varargin{i+1};
            case 'AL_ITERS'
                AL_iters = round(varargin{i+1});
                if (AL_iters <= 0 )
                    error('AL_iters must a positive integer');
                end
            case 'POSITIVITY'
                positivity = varargin{i+1};
                if strcmp(positivity,'yes')
                    reg_pos = 1;
                end
            case 'MU'
                mu = varargin{i+1};
                if mu <= 0
                    error('mu must be positive');
                end
            case 'VERBOSE'
                verbose = varargin{i+1};
            otherwise
                % Something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end


% test for image size correctness
if reg_TV > 0
    if N ~= prod(im_size)
        error('wrong image size')
    end
    n_lin = im_size(1);
    n_col = im_size(2);
    
    % build handlers and necessary stuff
    % horizontal difference operators
    FDh = zeros(im_size);
    FDh(1,1) = -1;
    FDh(1,end) = 1;
    FDh = fft2(FDh);
    FDhH = conj(FDh);
    
    % vertical difference operator
    FDv = zeros(im_size);
    FDv(1,1) = -1;
    FDv(end,1) = 1;
    FDv = fft2(FDv);
    FDvH = conj(FDv);
    
    IL = 1./( FDhH.* FDh + FDvH.* FDv + 1);
    
    Dh = @(x) real(ifft2(fft2(x).*FDh));
    DhH = @(x) real(ifft2(fft2(x).*FDhH));
    
    Dv = @(x) real(ifft2(fft2(x).*FDv));
    DvH = @(x) real(ifft2(fft2(x).*FDvH));
    
end


%%
%---------------------------------------------
% just least squares
%---------------------------------------------
if ~reg_TV && ~reg_wns && ~reg_pos
    X = pinv(A)*Y;
    res = norm(A*X-Y,'fro');
    return
end


%%
%---------------------------------------------
%  Constants and initializations
%---------------------------------------------

% number of regularizers
n_reg =  reg_wns + reg_pos + reg_TV;


%%
%---------------------------------------------
%  Initializations
%---------------------------------------------

% no intial solution supplied
if X0 == 0
    X = (A'*A + n_reg*eye(r))\A'*Y;
end

% what regularizers ?
%  1 - data term
%  2 - positivity
%  3 - weighted nucleus sparsity
%  4 - TV

index = 1

% initialize V variables
V = cell(1 + n_reg,1);

% initialize D variables (scaled Lagrange Multipliers)
D = cell(1 + n_reg,1);

%  data term (always present)
reg(1) = 1;             % regularizers
V{index} = A*X;         % V1
D{1} = zeros(size(Y));  % Lagrange multipliers

% next V
index = index + 1;
% POSITIVITY
if reg_pos == 1
    reg(index) = 2;
    V{index} = X;
    D{index} = zeros(size(X));
    index = index +1;
end
% Weighted nucleus sparsity
if reg_wns == 1
    reg(index) = 3;
    V{index} = X;
    D{index} = zeros(size(X));
    index = index +1;
end
% TV
if reg_TV == 1
    reg(index) = 4;
    V{index} = X;
    D{index} = zeros(size(X));
    
    % convert X into a cube
    X_im = reshape(X',im_size(1), im_size(2),r);
    
    % create two images per band (horizontal and vertical differences)
    V{index+1} = cell(r,2);
    D{index+1} = cell(r,2);
    for i=1:r
        V{index+1}{i}{1} = Dh(X_im(:,:,i));   % horizontal differences
        V{index+1}{i}{2} = Dv(X_im(:,:,i));   % horizontal differences
        D{index+1}{i}{1} = zeros(im_size);   % horizontal differences
        D{index+1}{i}{2} = zeros(im_size);   % horizontal differences
    end
    clear X_im;
end




%%
%---------------------------------------------
%  AL iterations - main body
%---------------------------------------------
tol1 = sqrt(N)*1e-6;
i=1;
res = inf;
while (i <= AL_iters) && (sum(abs(res)) > tol1)
    
    
    % solve the quadratic step (all terms depending on X)
    Xi = A'*(V{1}+D{1});
    for j = 2:(n_reg+1)
        Xi = Xi+ V{j} + D{j};
    end
    X = (A'*A + n_reg*eye(r))\Xi;

    % Weight vector of nucleus sparsity
    gw = [0*ones(1,im_size(1)*im_size(2));exp(-X(2,:));0*ones(1,im_size(1)*im_size(2));0*ones(1,im_size(1)*im_size(2))];

    % Compute the Mourau proximity operators
    for j=1:(n_reg+1)
        %  data term (V1)
        if  reg(j) == 1
            V{j} = (1/(1+mu)*(Y+mu*(A*X-D{j})));
        end
        %  positivity   (V2)
        if  reg(j) == 2
            V{j} = max(X-D{j},0);
        end
        % l1 norm  (V4)
        if  reg(j) == 3
            V{j} = soft(X-D{j}, (lambda_1/mu).*gw );
        end
        % TV  (V5 and V6)
        if  reg(j) == 4
            % update V5: solves the problem:
            %    min 0.5*||L*V5-(V6+D7)||^2+0.5*||V5-(X-d5)||^2
            %      V5
            %
            % update V6: min 0.5*||V6-(L*V5-D6)||^2 + lambda_tv * |||V6||_{1,1}
            
            nu_aux = X - D{j};
            % convert nu_aux into image planes
            % convert X into a cube
            nu_aux5_im = reshape(nu_aux',im_size(1), im_size(2),r);
            % compute V5 in the form of image planes
            for k =1:r
                % V5
                V5_im(:,:,k) = real(ifft2(IL.*fft2(DhH(V{j+1}{k}{1}+D{j+1}{k}{1}) ...
                    +  DvH(V{j+1}{k}{2}+D{j+1}{k}{2}) +  nu_aux5_im(:,:,k))));
                % V6
                aux_h = Dh(V5_im(:,:,k));
                aux_v = Dv(V5_im(:,:,k));
                V{j+1}{k}{1} = soft(aux_h - D{j+1}{k}{1}, (lambda_TV/mu) );   %horizontal
                V{j+1}{k}{2} = soft(aux_v - D{j+1}{k}{2}, (lambda_TV/mu) );   %vertical
                % update D6
                D{j+1}{k}{1} =  D{j+1}{k}{1} - (aux_h - V{j+1}{k}{1});
                D{j+1}{k}{2} =  D{j+1}{k}{2} - (aux_v - V{j+1}{k}{2});
            end
            % convert V6 to matrix format
            V{j} = reshape(V5_im, prod(im_size),r)';
            
        end
        
    end
    
    
    
    % update Lagrange multipliers
    
    for j=1:(n_reg+1)
        if  reg(j) == 1
            D{j} = D{j} - (A*X-V{j});
        else
            D{j} = D{j} - (X-V{j});
        end
    end
    
    % compute residuals
    for j=1:(n_reg+1)
        if  reg(j) == 1
            res(j) = norm(A*X-V{j},'fro');
        end
    end
    Res(end+1) = sum(abs(res));



    if mod(i,10) == 1
       
        % compute residuals
        st = [];
        for j=1:(n_reg+1)
            if  reg(j) == 1
                st = strcat(st,sprintf(' res(%i) = %2.6f',reg(j),res(j) ));
            else
                res(j) = norm(X-V{j},'fro');
                st = strcat(st,sprintf('  res(%i) = %2.6f',reg(j),res(j) ));
            end
        end
        if  strcmp(verbose,'yes')
            fprintf(strcat(sprintf('iter = %i - L1 = %d - LTV = %d -',i,lambda_1, lambda_TV ),st,'\n'));
        end

    end
    
    
    i=i+1;
    
end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %