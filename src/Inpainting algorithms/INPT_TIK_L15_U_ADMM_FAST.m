function [u,itrs] = INPT_TIK_L15_U_ADMM_FAST(b,M,mu,beta_r,u,itrs_th,itrs_rel_chg_th)
%
% Compute approximate solutions (that is, minimizers) of the
% unconstrained TIK-L1 variational model: 
%
% u*(mu) = argmin J(u;mu),   J(u;mu) = (1/2)||D u||_2^2 + mu ||S u - b||_1  
%             u
%
% for the non-blind inpainting, that is denoising + holes filling (non-blind 
% means that the selection/masking operator is assumed to be known), 
% of grayscale images using the ADMM (iterative) optimization algorithm

% INPUTS:
    % inputs defining the unconstrained TIK-L1 cost function to be minimized:
        % b         --> observed corrupted (masked and noisy) image
        % M         --> inpainting (binary) mask image
        % mu        --> regularization parameter (positive real scalar)        
    % inputs defining the ADMM (iterative) algorithm:
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
    d     = h * w;   % total number of pixels
    % d_M   = sum(M(:));

    % initialize optimization (primal and dual) variables
        % primal variables
            % the initial u is an input of the algorithm (of this Matlab function)
            % the initial r is not necessary (see below)!
            gamma = (1.5 * beta_r) / mu;
        % dual variables (Lagrange multipliers)
            % associated with the introduced auxiliary variable r
            lambda_r = zeros(h,w);

    % given the inputs  b, M, mu and beta_r,  for efficiency purposes 
    % compute/store once for all at the beginning (here) all the quantities 
    % that will be used for the solution of the ADMM algorithm sub-problems 
    % but that do not change during the ADMM iterations
    
        % sub-problem for the primal variable r (shrinkage)
        %mu_over_beta_r  = mu / beta_r; % constant in the shrinkage
        
        % sub-problem for the primal variable u (linear system)
          %
          I_DFT =  ones(h); 
  
          Dh_DFT     = psf2otf([1,-1],size(b));   Dv_DFT          = psf2otf([1;-1],size(b));

          DhT_DFT    = conj(Dh_DFT);              DvT_DFT         = conj(Dv_DFT);

          SUB_X_CM_DFT   = (1 ./ beta_r) * ( DhT_DFT .* Dh_DFT + DvT_DFT .* Dv_DFT ) + I_DFT;
        
        % sub-problem for the dual variable lambda_r
        % ... nothing to pre-compute!
            
% -------------------------
% carry out ADMM iterations
% -------------------------

% initialize iteration index and stopping criteria flags
itrs        = 0;   
stop_flags  = [0, 0];

% while stopping criteria are not satisfied...iterate

while ( ( sum(stop_flags) == 0 ) )
    
    % update iteration index
      itrs = itrs + 1;
    
    % ----------------------------------------------------------------
    % solve the ADMM sub-problem for the primal variable r (shrinkage)
    % ----------------------------------------------------------------
      %
      q    =  u + ( lambda_r ./ beta_r );

      q_1  =  b - q;
           
      r    = q - M .* sign(q_1) .* ( 9 / (8 * gamma^2) ) .* ( 1 - sqrt( 1 + ( 16*abs(q_1)*gamma^2 ) / 9 ) );
    
    % ---------------------------------------------------------------------
    % solve the ADMM sub-problem for the primal variable u (linear system)
    % ---------------------------------------------------------------------
        
        % store u_old for future computation of iterates relative change
          u_old     = u; 

        % compute the right-hand side of the linear system
         
          rhs   = ( r - lambda_r ./ beta_r );
        

        % compute the new iterate
          %
          u_DFT    = fft2(rhs) ./ SUB_X_CM_DFT;

          u        = real( ifft2( u_DFT ) );
          
         
                   
    % ---------------------------------------------------
    % compute u relative change (for stopping criterion)
    % ---------------------------------------------------
      %
      rel_chg = norm (u - u_old,'fro') / norm (u_old,'fro');
    
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
    
        % % gradient of the current iterate u
        %   [Dhu,Dvu] = D(u);
        % % residual image associated with the current iterate u
        %Su_b      = M .* u - b;
    
    % ------------------------------------------------------------------
    % update dual variables (Lagrange multipliers) lambdas (dual ascent)
    % ------------------------------------------------------------------
        
        % Lagrange multipliers associated with the auxiliary variable r
          %lambda_r = lambda_r - beta_r .* real( ifft2( r - u_DFT ) )

          lambda_r = lambda_r - beta_r .* ( r - u );
          
    % --------------------------------------
    % (eventually) compute/store/show debugs
    % --------------------------------------

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



