function [u,itrs] = REST_TV_L1_U_ADMM(b,blur_k,mu,beta_t,beta_r,u,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true)
%
% Compute approximate solutions (that is, minimizers) of the
% unconstrained TV-L1 variational model: 
%
% u*(mu) = argmin J(u;mu),   J(u;mu) = TV(u) + mu ||K u - b||_1  
%             u
%
% for the non-blind restoration, that is denoising + deblurring (non-blind 
% means that the blur operator is assumed to be known), of grayscale images 
% using the ADMM (iterative) optimization algorithm, as illustrated in [1].
%
% [1] Min Tao and Junfeng Yang, Alternating Direction Algorithms for 
%     Total Variation Deconvolution in Image Reconstruction, TR0918, 
%     Department of Mathematics, Nanjing University, 2009.

% INPUTS:
    % inputs defining the unconstrained TV-L1 cost function to be minimized:
        % b         --> observed corrupted (blurred and noisy) image
        % blur_k    --> blur PSF (kernel of spatial 2D convolution)
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
%     to restore, that is both the blur and the first-order derivatives 
%     near the image boundaries are computed according to this assumption!
%     This assumption allows to diagonalize the blur matrix K and the 
%     first-order horizontal and vertical partial derivative matrices
%     Dh and Dv by means of the 2D Discrete Fourier Transform (DFT)
    
% --------------------
% initialize algorithm
% --------------------
    
    % extract image dimensions
    [h,w] = size(b); % height and width, in pixels
    d     = h * w;   % total number of pixels
    
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

    % given the inputs  b, blur_k, mu, beta_t and beta_r,  for efficiency purposes 
    % compute/store once for all at the beginning (here) all the quantities 
    % that will be used for the solution of the ADMM algorithm sub-problems 
    % but that do not change during the ADMM iterations
    
        % sub-problem for the primal variable t (shrinkage)
        one_over_beta_t = 1 / beta_t; % constant in the shrinkage
        
        % sub-problem for the primal variable r (shrinkage)
        mu_over_beta_r  = mu / beta_r; % constant in the shrinkage
        
        % sub-problem for the primal variable u (linear system)
        beta_r_over_beta_t = beta_r / beta_t; % constant in the coefficient matrix and in the right-hand side
        K_DFT              = psf2otf(blur_k,[h,w]); % 2D DFT for blur matrix K
        Dh_DFT             = psf2otf([1,-1],[h,w]); % 2D DFT for finite difference matrix Dh
        Dv_DFT             = psf2otf([1;-1],[h,w]); % 2D DFT for finite difference matrix Dv
        KT_DFT             = conj(K_DFT);  % 2D DFT for transpose of blur matrix K
        DhT_DFT            = conj(Dh_DFT); % 2D DFT for transpose of finite difference matrix Dh
        DvT_DFT            = conj(Dv_DFT); % 2D DFT for transpose of finite difference matrix Dv
        SUB_X_CM_DFT       = DhT_DFT .* Dh_DFT + DvT_DFT .* Dv_DFT + beta_r_over_beta_t * (KT_DFT .* K_DFT);
        beta_r_over_beta_t_KT_DFT = beta_r_over_beta_t * KT_DFT;
        
        % sub-problem for the dual variables lambda_t and lambda_r
        % ... nothing to pre-compute!
            
    % initialize other necessary quantities
        % gradient of the initial iterate u
        [Dhu,Dvu] = D(u); 
        % residual image associated with the initial iterate u
        Ku_b      = real(ifft2( K_DFT .* fft2(u) )) - b;
    
% --------------------------------------
% (eventually) compute/store/show debugs
% --------------------------------------
if ( itrs_debugs == 1 )
    % allocate arrays where storing iterations-based debug quantities
    J_regs        = zeros(1,itrs_th+1);
    J_fids        = zeros(1,itrs_th+1);
    rel_chgs      = zeros(1,itrs_th);
    Ku_b_means    = zeros(1,itrs_th+1);
    taus          = zeros(1,itrs_th+1);
    snrs          = zeros(1,itrs_th+1);
    ssims         = zeros(1,itrs_th+1);
    % compute/initialize and store iterations-based debug quantities
    J_regs(1)     = sum( sqrt(Dhu(:).^2 + Dvu(:).^2) );
    J_fids(1)     = mu * sum( abs(Ku_b(:)) );
    Ku_b_means(1) = mean( Ku_b(:) );
    taus(1)       = sqrt( sum( Ku_b(:).^2 ) ) / (sqrt(d)*n_sigma);
    snrs(1)       = compute_snr(u,u_true);
    ssims(1)      = ssim(u,u_true);  
    % Matlab command window debugs
    fprintf('\n---------------------------------------------------------');
    fprintf('\nTV-L1-U-ADMM (mu = %7.2f) ITERATIONS-BASED DEBUGS:',mu);
    fprintf('\n---------------------------------------------------------');
    fprintf('\nIT %04d:  RC = %9.7f, J = %13.6f, tau* = %6.4f, [ISNR,ISSIM] = [%6.3f , %6.4f]',...
        0,0,(J_regs(1)+J_fids(1)),taus(1),snrs(1)-snrs(1),ssims(1)-ssims(1));
end
    

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
    % ----------------------------------------------------------------
    q_h                 = Dhu + lambda_t_h;
    q_v                 = Dvu + lambda_t_v;
    q_norm              = sqrt( q_h.^2 + q_v.^2 );
    q_norm(q_norm == 0) = one_over_beta_t;
    q_norm              = max( q_norm - one_over_beta_t , 0 ) ./ q_norm;
    t_h                 = q_norm .* q_h;
    t_v                 = q_norm .* q_v;
    
    % ----------------------------------------------------------------
    % solve the ADMM sub-problem for the primal variable r (shrinkage)
    % ----------------------------------------------------------------
    w = Ku_b + lambda_r;
    r = sign(w) .* max( abs(w) - mu_over_beta_r , 0 );
    
    % --------------------------------------------------------------------
    % solve the ADMM sub-problem for the primal variable u (linear system)
    % --------------------------------------------------------------------
        
        % store u_old for future computation of iterates relative change
        u_old       = u; 
        % compute the (DFT of the) right-hand side of the linear system
        rhs_DFT     = fft2( DT(t_h - lambda_t_h,t_v - lambda_t_v) ) + ... 
                      beta_r_over_beta_t_KT_DFT .* fft2(r - lambda_r + b);
        % compute the (DFT of the) new iterate
        u_DFT       = rhs_DFT ./ SUB_X_CM_DFT;
        % compute the new iterate
        u           = real(ifft2(u_DFT));
                   
    % --------------------------------------------------
    % compute u relative change (for stopping criterion)
    % --------------------------------------------------
    rel_chg = norm(u - u_old,'fro') / norm(u_old,'fro');
    
    % -----------------------
    % check stopping criteria
    % -----------------------
    if ( itrs == itrs_th )
        stop_flags(1) = 1;
    end
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
        Ku_b      = real(ifft2( K_DFT .* u_DFT )) - b;
    
    % ------------------------------------------------------------------
    % update dual variables (Lagrange multipliers) lambdas (dual ascent)
    % ------------------------------------------------------------------
        
        % Lagrange multipliers associated with the auxiliary variable t = (th;tv)
        lambda_t_h = lambda_t_h - (t_h - Dhu);
        lambda_t_v = lambda_t_v - (t_v - Dvu);
        % Lagrange multipliers associated with the auxiliary variable r
        lambda_r   = lambda_r   - (r - Ku_b);
        
    % --------------------------------------
    % (eventually) compute/store/show debugs
    % --------------------------------------
    if ( itrs_debugs == 1 )
        % compute and store iterations-based debug quantities
        J_regs(itrs+1)     = sum( sqrt(Dhu(:).^2 + Dvu(:).^2) );
        J_fids(itrs+1)     = mu * sum( abs(Ku_b(:)) );
        rel_chgs(itrs)     = rel_chg;
        Ku_b_means(itrs+1) = mean( Ku_b(:) );
        taus(itrs+1)       = sqrt( sum( Ku_b(:).^2 ) ) / (sqrt(d)*n_sigma);
        snrs(itrs+1)       = compute_snr(u,u_true);
        ssims(itrs+1)      = ssim(u,u_true);
        % Matlab command window debugs
        fprintf('\nIT %04d:  RC = %9.7f, J = %13.6f, tau* = %6.4f, [ISNR,ISSIM] = [%6.3f , %6.4f]',...
            itrs,rel_chgs(itrs),(J_regs(itrs+1)+J_fids(itrs+1)),taus(itrs+1),snrs(itrs+1)-snrs(1),ssims(itrs+1)-ssims(1));
    end
    
end

% one of the two stopping criteria satisfied --> end iterations  
% and, before returning the restored image u and the total number 
% of performed ADMM iterations itrs to the main script, show 
% (eventually if itrs_debugs = 1) some iterations-based debug quantities

% show iterations-based debug quantities
if ( itrs_debugs == 1 )
    
    % first, shrink arrays where we stored iterations-based debug quantities
    J_regs        = J_regs(1:itrs+1);
    J_fids        = J_fids(1:itrs+1);
    rel_chgs      = rel_chgs(1:itrs);
    Ku_b_means    = Ku_b_means(1:itrs+1);
    taus          = taus(1:itrs+1);
    snrs          = snrs(1:itrs+1);
    ssims         = ssims(1:itrs+1);
    
    % then show all debugs
    fprintf('\n');
    
    plots_x_min = 0;
    plots_x_max = itrs;
    
    figure('Name','Figure 40:  TV-L1-U-ADMM:   ITERATIONS-BASED DEBUGS','NumberTitle','off');
    set(gcf,'Position',get(0,'ScreenSize'));
    
    subplot(2,3,1)
    plot(0:itrs,J_regs,'b');
    xlabel('k  (ADMM iteration index)');
    ylabel('$\mathrm{TV}\left(u^{(k)}\right) = \sum_i \left\| (\nabla u)_i \right\|_2$','interpreter','latex');
    plot1_y_min = min( J_regs(2:end) );
    plot1_y_max = max( J_regs(2:end) );
    plot1_y_min = plot1_y_min - 0.05 * (plot1_y_max - plot1_y_min);
    plot1_y_max = plot1_y_max + 0.05 * (plot1_y_max - plot1_y_min);
    axis([plots_x_min,plots_x_max,plot1_y_min,plot1_y_max]);
    
    subplot(2,3,2)
    plot(0:itrs,J_fids,'k');
    xlabel('k  (ADMM iteration index)');
    ylabel('$\mu \,\, \mathrm{L}_1\left(u^{(k)};b,K\right) = \mu \,\, \left\|Ku^{(k)}-b\right\|_1$','interpreter','latex');
    plot2_y_min = min( J_fids(2:end) );
    plot2_y_max = max( J_fids(2:end) );
    plot2_y_min = plot2_y_min - 0.05 * (plot2_y_max - plot2_y_min);
    plot2_y_max = plot2_y_max + 0.05 * (plot2_y_max - plot2_y_min);
    axis([plots_x_min,plots_x_max,plot2_y_min,plot2_y_max]);
    
    subplot(2,3,3)
    plot(0:itrs,J_regs + J_fids,'r');
    xlabel('k  (ADMM iteration index)');
    ylabel('$J\left(u^{(k)};\mu\right) = \mathrm{TV}\left(u^{(k)}\right) + \mu \, \mathrm{L}_1\left(u^{(k)};b,K\right)$','interpreter','latex');
    plot3_y_min = min( J_regs(2:end) + J_fids(2:end) );
    plot3_y_max = max( J_regs(2:end) + J_fids(2:end) );
    plot3_y_min = plot3_y_min - 0.05 * (plot3_y_max - plot3_y_min);
    plot3_y_max = plot3_y_max + 0.05 * (plot3_y_max - plot3_y_min);
    axis([plots_x_min,plots_x_max,plot3_y_min,plot3_y_max]);
    
    subplot(2,3,4)
    plot(0:(itrs-1),rel_chgs,'b');
    xlabel('k  (ADMM iteration index)');
    ylabel('$\delta^{(k)} = \left\| u^{(k+1)}-u^{(k)} \right\|_2 / \left\|u^{(k)}\right\|_2$','interpreter','latex');
    plot4_y_min = min( rel_chgs(2:end) );
    plot4_y_max = max( rel_chgs(2:end) );
    plot4_y_min = plot4_y_min - 0.05 * (plot4_y_max - plot4_y_min);
    plot4_y_max = plot4_y_max + 0.05 * (plot4_y_max - plot4_y_min);
    axis([plots_x_min,plots_x_max,plot4_y_min,plot4_y_max]);
    
    subplot(2,3,5)
    plot(0:itrs,taus,'b');
    xlabel('k  (ADMM iteration index)');
    ylabel('$\tau^{(k)} = \left\|K u^{(k)} - b\right\|_2 / (\sqrt{d} \: \sigma_n)$','interpreter','latex');
    plot5_y_min = min( taus(2:end) );
    plot5_y_max = max( taus(2:end) );
    plot5_y_min = plot5_y_min - 0.05 * (plot5_y_max - plot5_y_min);
    plot5_y_max = plot5_y_max + 0.05 * (plot5_y_max - plot5_y_min);
    axis([plots_x_min,plots_x_max,plot5_y_min,plot5_y_max]);
    
    subplot(2,3,6)    
    yyaxis left;
    plot(0:itrs,snrs-snrs(1));
    xlabel('k  (ADMM iteration index)'); 
    ylabel('$ISNR\left(u_{true},u^{(k)}\right)$','interpreter','latex');
    plot6_y_min = min( snrs(2:end)-snrs(1) );
    plot6_y_max = max( snrs(2:end)-snrs(1) );
    plot6_y_max = plot6_y_max + 0.05 * (plot6_y_max - plot6_y_min);
    axis([plots_x_min,plots_x_max,plot6_y_min,plot6_y_max]);
    yyaxis right;
    plot(0:itrs,ssims-ssims(1));
    xlabel('k  (ADMM iteration index)');
    ylabel('$ISSIM\left(u_{true},u^{(k)}\right)$','interpreter','latex');
    plot6_y_min = min( ssims(2:end)-ssims(1) );
    plot6_y_max = max( ssims(2:end)-ssims(1) );
    plot6_y_max = plot6_y_max + 0.05 * (plot6_y_max - plot6_y_min);
    axis([plots_x_min,plots_x_max,plot6_y_min,plot6_y_max]);
    
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


