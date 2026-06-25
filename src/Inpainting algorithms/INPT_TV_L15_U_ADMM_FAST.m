function [u,itrs] = INPT_TV_L15_U_ADMM_FAST(b,M,mu,beta_t,beta_r,u,itrs_th,itrs_rel_chg_th)
%
% Compute approximate solutions (that is, minimizers) of the
% unconstrained TV-L1 variational model: 
%
% u*(mu) = argmin J(u;mu),   J(u;mu) = TV(u) + mu ||S u - b||_1  
%             u
%
% for the non-blind inpainting, that is denoising + holes filling (non-blind 
% means that the selection/masking operator is assumed to be known), 
% of grayscale images using the ADMM (iterative) optimization algorithm

% INPUTS:
    % inputs defining the unconstrained TV-L1 cost function to be minimized:
        % b         --> observed corrupted (masked and noisy) image
        % M         --> inpainting (binary) mask image
        % mu        --> regularization parameter (positive real scalar)        
    % inputs defining the ADMM (iterative) algorithm:
        % beta_t    --> penalty parameter associated with the auxiliary variable t (positive real scalar)
        % beta_r    --> penalty parameter associated with the auxiliary variable r (positive real scalar)
        % u         --> initial iterate (usually, but not necessarily, the observed corrupted image b)
        % ADMM stopping criteria:
            % itrs_th         --> threshold: maximum number of ADMM iterations
            % itrs_rel_chg_th --> threshold: iterates relative change
    % (optional) inputs defining eventual "debugs" 
        % itrs_debugs --> flag: 0->no debug; 1->debug
        % n_sigma     --> noise standard deviation
        % u_true      --> original uncorrupted image
        
% We remark the following:
%   - all images are stored as matrices, that is they are not vectorized!
%   - periodic boundary conditions are assumed for the unknown clean image 
%     to restore, that is the first-order derivatives near the image 
%     boundaries are computed according to this assumption!
    
% --------------------
% initialize algorithm
% --------------------
    
    % extract image dimensions
    [h,w] = size(b); % height and width, in pixels
    %d     = h * w;   % total number of pixels
    %d_M   = sum(M(:));
    
    % initialize optimization (primal and dual) variables
        % primal variables
            % the initial u is an input of the algorithm (of this Matlab function)
            % the initial t is not necessary (see below)!
            % the initial r is not necessary (see below)!

        % dual variables (Lagrange multipliers)
            % associated with the introduced auxiliary variable t = (th;tv)
              lambda_t_h = zeros(h,w);
              lambda_t_v = zeros(h,w);

            % associated with the introduced auxiliary variable r
              lambda_r   = zeros(h,w);
     

    % given the inputs  b, M, mu, beta_t and beta_r,  for efficiency purposes 
    % compute/store once for all at the beginning (here) all the quantities 
    % that will be used for the solution of the ADMM algorithm sub-problems 
    % but that do not change during the ADMM iterations
    
        % sub-problem for the primal variable t (shrinkage)
          one_over_beta_t = 1 / beta_t; % constant in the shrinkage
        
        % sub-problem for the primal variable r (shrinkage)
          gamma = (1.5 * beta_r) / mu;
        
        
        % sub-problem for the primal variable u (linear system)
          beta_r_over_beta_t = beta_r / beta_t; % constant in the coefficient matrix and in the right-hand side

        I_DFT  = ones(h);

        Dh_DFT  = psf2otf( [1, -1], size(b) );  Dv_DFT = psf2otf( [1; -1], size(b) );

        DhT_DFT = conj( Dh_DFT );              DvT_DFT = conj( Dv_DFT );

        SUB_X_CM_DFT     = ( DhT_DFT .* Dh_DFT + DvT_DFT .* Dv_DFT )  + beta_r_over_beta_t .* I_DFT;
        
        % sub-problem for the dual variables lambda_t and lambda_r
        % ... nothing to pre-compute!
            
    % initialize other necessary quantities
        % gradient of the initial iterate u
        % [Dhu,Dvu] = D(u); 
        % residual image associated with the initial iterate u
        % Su_b      = M .* u - b;

% initialize other necessary quantities
% gradient of the initial iterate u
[Dhu,Dvu] = D(u); 

% -------------------------
% carry out ADMM iterations
% -------------------------

% initialize iteration index and stopping criteria flags
itrs        = 0;   
stop_flags  = [0,0];

% while stopping criteria are not satisfied...iterate
while ( sum(stop_flags) == 0 )
    
    % update iteration index
    itrs = itrs + 1;
    
    % ----------------------------------------------------------------
    % solve the ADMM sub-problem for the primal variable t (shrinkage)
    % ---------------------------------------------------------------

    q_h                 = Dhu + lambda_t_h ./ beta_t;
    q_v                 = Dvu + lambda_t_v ./ beta_t;

    q_norm              = sqrt( q_h.^2 + q_v.^2 );
    q_norm(q_norm == 0) = one_over_beta_t;
    q_norm              = max( q_norm - one_over_beta_t , 0 ) ./ q_norm;

    t_h                 = q_norm .* q_h;

    t_v                 = q_norm .* q_v;
    
    % ----------------------------------------------------------------
    % solve the ADMM sub-problem for the primal variable r (shrinkage)
    % ---------------------------------------------------------------
      %
      q   = u + (lambda_r ./ beta_r );

      q_1 = b - q;

      r   = q - M .* ( 9 / (8 * gamma^2) ) .* sign(q_1) .* ( 1 - sqrt( 1 + ( 16*abs(q_1)*gamma^2 ) / 9 ) );
    
    % --------------------------------------------------------------------
    % solve the ADMM sub-problem for the primal variable u (linear system)
    % --------------------------------------------------------------------
        
      % store u_old for future computation of iterates relative change
        u_old = u;  

      % compute the right-hand side of the linear system
        %    

        rhs  =  DT( t_h - lambda_t_h ./ beta_t, t_v - lambda_t_v ./ beta_t ) + ...
                    (beta_r_over_beta_t) .* (r - lambda_r ./ beta_r);

        rhs_DFT = fft2( rhs );

       
        % compute the new iterate
          %
          u_DFT = rhs_DFT ./ SUB_X_CM_DFT;

          u = real( ifft2 ( u_DFT ) );
                   
    % --------------------------------------------------
    % compute u relative change (for stopping criterion)
    % --------------------------------------------------
      rel_chg = norm(u - u_old, 'fro') / norm(u_old, 'fro');
    
    % -----------------------
    % check stopping criteria
    % -----------------------
    if ( itrs == itrs_th )
        stop_flags(1) = 1;

    end
    %
    if ( rel_chg < itrs_rel_chg_th )  
        stop_flags(2) = 1;

    end
    
    % -------------------------------
    % compute quantities used for the 
    % subsequent lambdas update step
    % and for the next iteration
    % -------------------------------
    
        % gradient of the current iterate u
        [Dhu,Dvu] = D(u);
        % residual image associated with the current iterate u
        % Su_b      = M .* u - b;
    
    % ------------------------------------------------------------------
    % update dual variables (Lagrange multipliers) lambdas (dual ascent)
    % ------------------------------------------------------------------
        
        % Lagrange multipliers associated with the auxiliary variable t = (th;tv)
          lambda_t_h   = lambda_t_h - beta_t .* (t_h - Dhu );
          lambda_t_v   = lambda_t_v - beta_t .* (t_v - Dvu );

        % Lagrange multipliers associated with the auxiliary variable r
          lambda_r   = lambda_r - beta_r .* (r - u);
   
end

%-------------------------
% NESTED FUNCTIONS
%-------------------------
function [Dhu,Dvu] = D(u)
%
% Given the image u, compute Du , that is the gradient of u:
%       Du = (Dh u ; Dv u) ,
% with horizontal and vertical partial derivates (Dh u) and (Dv u) 
% both discretized by forward finite differences (kernel [-1,1])

Dhu = [ u(:,2:end) - u(:,1:(end-1)) , u(:,1) - u(:,end) ];
Dvu = [ u(2:end,:) - u(1:(end-1),:) ; u(1,:) - u(end,:) ]; 

end

function [DT_u12] = DT(u1,u2)
%
% Given the two images u1, u2, compute DT(u1;u2), that is:
%       DT(u1;u2) = Dh^T u1 + Dv^T u2,
% with horizontal and vertical partial derivates (Dh u) and (Dv u) 
% both discretized by forward finite differences (kernel [-1,1])

DhT_u1 = [ u1(:,end) - u1(:,1) , -diff(u1,1,2) ];
DvT_u2 = [ u2(end,:) - u2(1,:) ; -diff(u2,1,1) ];
DT_u12 = DhT_u1 + DvT_u2;

end


end


