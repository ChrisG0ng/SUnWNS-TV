function [X, Res] = SUnWNS_TV_opt(A, Y, varargin)
% SUnWNS_TV (optimized) — faster version 
%  * Same signature/outputs, improved speed & stability on CPU
%
% Key changes:
%  - Switch between weighted nuclues sparsity and L1 sparsity
%  - Pre‑factorize (A'*A + n_reg*I) with Cholesky and reuse each iter
%  - Fix WNS weight size bug (works even when TV is off)
%  - TV subproblem builds RHS fully in Fourier domain (fewer fft/ifft)
%  - Preallocate large arrays; replace dynamic growth of Res
%  - ADMM stop (primal/dual residuals) instead of fixed tol on data term
%  - Optional knobs: PRECISION (double/single), FFTW_PLANNER ('estimate'|'measure'|'patient')
%
% Usage is unchanged; optional new params (safe defaults):
%  'PRECISION'    {'double','single'}  default 'single'
%  'FFTW_PLANNER' {'estimate','measure','patient'} default 'estimate'
%% ---------------- argument checks ----------------
if (nargin-length(varargin)) ~= 2
    error('Wrong number of required parameters');
end
[LM,r] = size(A); [L,N] = size(Y);
if (LM ~= L), error('mixing matrix A and data set Y are inconsistent'); end
%% ---------------- defaults ----------------
sparsity_mode = 'wns';   % 'wns' (exp(-Xh)) or 'l1' (gw = 1)
reg_wns  = 0;        % L1 on H row
reg_TV   = 0;        % TV
im_size  = [];       % [h, w]
AL_iters = 1000;     % max ADMM iters
mu       = 0.005;     % AL weight
verbose  = 'off';    % 'yes' to print
reg_pos  = 0;        % positivity
X0       = 0;        % no warm start
Res      = [];       % residual history (will be preallocated later)
% New (optional) tuning knobs — backwards compatible
precision = 'single';    % or 'double'
fftw_planner = 'estimate';
%% ---------------- read optional parameters ----------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for iopt=1:2:(length(varargin)-1)
        switch upper(varargin{iopt})
            case 'SPARSITY_MODE'
                sparsity_mode = lower(varargin{iopt+1});
                if ~ismember(sparsity_mode, {'wns','l1'})
                    error('SPARSITY must be ''wns'' or ''l1''.');
                end
            case 'LAMBDA_1'
                lambda_1 = varargin{iopt+1};
                if lambda_1 < 0, error('lambda must be positive'); end
                if lambda_1 > 0, reg_wns = 1; end
            case 'LAMBDA_TV'
                lambda_TV = varargin{iopt+1};
                if lambda_TV < 0, error('lambda must be non-negative'); end
                if lambda_TV > 0, reg_TV = 1; end
            case 'IM_SIZE'
                im_size = varargin{iopt+1};
            case 'AL_ITERS'
                AL_iters = round(varargin{iopt+1});
                if (AL_iters <= 0 ), error('AL_iters must be positive'); end
            case 'POSITIVITY'
                positivity = varargin{iopt+1};
                if strcmp(positivity,'yes'), reg_pos = 1; end
            case 'MU'
                mu = varargin{iopt+1};
                if mu <= 0, error('mu must be positive'); end
            case 'VERBOSE'
                verbose = varargin{iopt+1};
            % --- new optional (safe if omitted) ---
            case 'PRECISION'
                precision = varargin{iopt+1};
            case 'FFTW_PLANNER'
                fftw_planner = varargin{iopt+1};
            otherwise
                error(['Unrecognized option: ''', varargin{iopt}, '''']);
        end
    end
end
% TV requires a consistent image size
if reg_TV
    if N ~= prod(im_size), error('wrong image size'); end
    n_lin = im_size(1); n_col = im_size(2);
end
% Cast precision once (CPU path)
if strcmpi(precision,'single')
    A = single(A); Y = single(Y); mu = single(mu);
    if reg_wns, lambda_1 = single(lambda_1); end
    if reg_TV,  lambda_TV = single(lambda_TV); end
else
    A = double(A); Y = double(Y); mu = double(mu);
    if reg_wns, lambda_1 = double(lambda_1); end
    if reg_TV,  lambda_TV = double(lambda_TV); end
end
% Optional FFTW plan tuning (CPU FFT backend)
try
    prevPlan = fftw('planner'); %#ok<NASGU>
    switch lower(fftw_planner)
        case {'estimate','measure','patient'}
            fftw('planner', lower(fftw_planner));
        otherwise
            % ignore invalid value, keep MATLAB default
    end
catch
    % older MATLAB: ignore
end
%% ---------------- prepare TV operators (as in your original) ----------------
if reg_TV > 0
    % horizontal diff frequency response
    FDh  = zeros(im_size, class(Y)); FDh(1,1) = -1; FDh(1,end) = 1; FDh  = fft2(FDh);
    FDhH = conj(FDh);
    % vertical diff frequency response
    FDv  = zeros(im_size, class(Y)); FDv(1,1) = -1; FDv(end,1) = 1; FDv  = fft2(FDv);
    FDvH = conj(FDv);
    % inverse filter for (Dh^H Dh + Dv^H Dv + I)
    IL   = 1 ./ ( FDhH.*FDh + FDvH.*FDv + 1 );
end
%% ---------------- LS warm-start & constants ----------------
n_reg = reg_wns + reg_pos + reg_TV;      % number of regularizers
if X0 == 0
    % LS with mild ridge (keeps SPD): use Cholesky and reuse per iter
    M = (A'*A + n_reg*eye(r, class(A)));
    [R,flag] = chol(M,'upper');
    if flag==0
        solveM = @(B) (R \ (R' \ B));
    else
        % robust fallback
        if exist('decomposition','builtin')
            dM = decomposition(M);
            solveM = @(B) (dM \ B);
        else
            solveM = @(B) (M \ B);
        end
    end
    X = solveM(A'*Y);
else
    X = X0;
end
% regularizer bookkeeping (preserve your original logic)
index = 1; 
reg = zeros(1 + n_reg,1);
V = cell(1 + n_reg,1);
D = cell(1 + n_reg,1);
% data term
reg(1) = 1; V{1} = A*X; D{1} = zeros(size(Y), class(Y));
idx = 2;
if reg_pos
    reg(idx) = 2; V{idx} = X; D{idx} = zeros(size(X), class(X)); idx = idx+1;
end
if reg_wns
    reg(idx) = 3; V{idx} = X; D{idx} = zeros(size(X), class(X)); idx = idx+1;
end
if reg_TV
    reg(idx) = 4; V{idx} = X; D{idx} = zeros(size(X), class(X));
    % convert X into cube & preallocate per band arrays (your layout)
    X_im = reshape(X', n_lin, n_col, r);
    V{idx+1} = cell(r,2); D{idx+1} = cell(r,2);
    % preallocate and initialize diffs/duals
    for k=1:r
        V{idx+1}{k}{1} = real(ifft2(fft2(X_im(:,:,k)).*FDh)); % Dh(X)
        V{idx+1}{k}{2} = real(ifft2(fft2(X_im(:,:,k)).*FDv)); % Dv(X)
        D{idx+1}{k}{1} = zeros(im_size, class(Y));
        D{idx+1}{k}{2} = zeros(im_size, class(Y));
    end
    clear X_im;
end
%% ---------------- ADMM main loop ----------------
Res = zeros(AL_iters,1, class(Y));   % preallocate history
% ADMM tolerances (Boyd)
tol_abs = sqrt(numel(X))*cast(1e-4,'like',Y);
tol_rel = cast(1e-3,'like',Y);
it = 1; resP = inf; resD = inf; 
Xprev = X;  % for dual residual
while it <= AL_iters
    % ---- X step ----
    Xi = A'*(V{1} + D{1});
    for j = 2:(n_reg+1)
        Xi = Xi + (V{j} + D{j});
    end
    X = solveM(Xi);
    % ---- WNS weights/L1 sparsity ----
    if reg_wns
        gw = zeros(r, N, 'like', X);
        if strcmp(sparsity_mode, 'wns')
            % weighted nucleus sparsity (default)
            gw(2,:) = exp(-X(2,:)); % exp(-xh): small xh -> larger weight
        else
            % plain L1 sparsity
            gw = 1;
        end
    end
    % ---- Prox steps ----
    for j = 1:(n_reg+1)
        switch reg(j)
            case 1  % data term (V1)
                V{j} = (Y + mu*(A*X - D{j})) / (1 + mu);
            case 2  % positivity (V?)
                V{j} = max(X - D{j}, 0);
            case 3  % weighted l1 on H row only
                V{j} = soft(X - D{j}, (lambda_1/mu).*gw );
            case 4  % TV (V5 and V6 in your comments)
                % nu_aux = X - D5
                nu_aux = X - D{j};
                nu_im  = reshape(nu_aux', n_lin, n_col, r);
                V5_im  = zeros(n_lin, n_col, r, 'like', Y); % preallocate once per iter
                for k=1:r
                    % Build RHS fully in frequency domain: fewer fft/ifft
                    RHS_hat = FDhH .* fft2(V{j+1}{k}{1} + D{j+1}{k}{1}) + ...
                              FDvH .* fft2(V{j+1}{k}{2} + D{j+1}{k}{2}) + ...
                              fft2(nu_im(:,:,k));
                    V5_im(:,:,k) = real(ifft2( IL .* RHS_hat ));
                    % Gradients of V5
                    V5_hat2 = fft2(V5_im(:,:,k));
                    aux_h = real(ifft2(FDh .* V5_hat2));
                    aux_v = real(ifft2(FDv .* V5_hat2));
                    % Prox for anisotropic TV (soft threshold)
                    V{j+1}{k}{1} = soft(aux_h - D{j+1}{k}{1}, (lambda_TV/mu));
                    V{j+1}{k}{2} = soft(aux_v - D{j+1}{k}{2}, (lambda_TV/mu));
                    % Dual updates for TV diffs
                    D{j+1}{k}{1} = D{j+1}{k}{1} - (aux_h - V{j+1}{k}{1});
                    D{j+1}{k}{2} = D{j+1}{k}{2} - (aux_v - V{j+1}{k}{2});
                end
                % Write back V5 in matrix form
                V{j} = reshape(V5_im, N, r)';
        end
    end
    % ---- Dual updates (all constraints) ----
    for j = 1:(n_reg+1)
        if reg(j) == 1
            D{j} = D{j} - (A*X - V{j});
        else
            D{j} = D{j} - (X - V{j});
        end
    end
    % ---- Residuals: primal & dual (Boyd) ----
    r_data = norm(A*X - V{1}, 'fro');
    r_pos = 0; r_wns = 0; r_tv = 0;
    if reg_pos
        idxp = find(reg==2,1);
        r_pos = norm(X - V{idxp}, 'fro');
    end
    if reg_wns
        idxw = find(reg==3,1);
        r_wns = norm(X - V{idxw}, 'fro');
    end
    if reg_TV
        idxt = find(reg==4,1);
        r_tv = norm(X - V{idxt}, 'fro');
    end
    resP = r_data + r_pos + r_wns + r_tv;
    resD = mu * norm(X - Xprev, 'fro');
    Res(it) = resP;
    % Stop rule
    eps_pri = tol_abs + tol_rel * max([norm(A*X,'fro'), norm(V{1},'fro')]);
    eps_dua = tol_abs + tol_rel * norm(mu*D{1},'fro');
    if (resP < 1*eps_pri && resD < 1*eps_dua)
        break;
    end
    Xprev = X;
    if strcmp(verbose,'yes') && mod(it,50)==1
        fprintf('it=%4d | r=%.3e s=%.3e\n', it, double(resP), double(resD));
    end
    it = it + 1;
end
% trim Res to actual length
Res = double(Res(1:it-1));
end  % ---- main function ----
% ---- local soft-threshold prox ----
function Z = soft(X, T)
% element-wise soft-threshold; T can be scalar or array (broadcasted)
Z = sign(X) .* max(abs(X) - T, 0);
end
