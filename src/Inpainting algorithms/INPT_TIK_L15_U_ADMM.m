function [u,itrs] = INPT_TIK_L15_U_ADMM(b,M,mu,beta_r,u,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true)
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
    d_M   = sum(M(:));
    
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
        dh           = ones(w,1);
        Dh_1d        = spdiags([-dh,dh],[0,1],w,w);
        Dh_1d(end,1) = 1;
        Dh           = kron( Dh_1d , speye(h) );

        dv           = ones(h,1);
        Dv_1d        = spdiags([-dv,dv],[0,1],h,h);
        Dv_1d(end,1) = 1;
        Dv           = kron( speye(w) , Dv_1d );

        SUB_X_CM     = (1 / beta_r) * ( Dh' * Dh + Dv' * Dv) + spdiags(M(:), 0, d, d);
        
        % sub-problem for the dual variable lambda_r
        % ... nothing to pre-compute!
            
    % initialize other necessary quantities
        % gradient of the initial iterate u
        [Dhu,Dvu] = D(u); 
        % residual image associated with the initial iterate u
        Su_b      = M .* u - b;
    
% --------------------------------------
% (eventually) compute/store/show debugs
% --------------------------------------
if ( itrs_debugs == 1 )
    % allocate arrays where storing iterations-based debug quantities
    J_regs        = zeros(1,itrs_th+1);
    J_fids        = zeros(1,itrs_th+1);
    rel_chgs      = zeros(1,itrs_th);
    Su_b_means    = zeros(1,itrs_th+1);
    taus          = zeros(1,itrs_th+1);
    snrs          = zeros(1,itrs_th+1);
    ssims         = zeros(1,itrs_th+1);
    % compute/initialize and store iterations-based debug quantities
    Su_b_vec      = Su_b( M > 0 );
    J_regs(1)     = (1/2) * sum( Dhu(:).^2 + Dvu(:).^2 );
    J_fids(1)     = mu * sum( abs(Su_b_vec) );
    Su_b_means(1) = mean( Su_b_vec );
    taus(1)       = sqrt( sum( Su_b_vec.^2 ) ) / (sqrt(d_M)*n_sigma);
    snrs(1)       = compute_snr(u,u_true);
    ssims(1)      = ssim(u,u_true); 
    % Matlab command window debugs
    fprintf('\n---------------------------------------------------------');
    fprintf('\nTIK-L1-U-ADMM (mu = %7.2f) ITERATIONS-BASED DEBUGS:',mu);
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
    % solve the ADMM sub-problem for the primal variable r (shrinkage)
    % ----------------------------------------------------------------
    q = Su_b + (1 / beta_r) .* lambda_r;

    r = q + 9 / (8*gamma^2) .* sign(q) .* ( 1 - sqrt(1 + (16*abs(q)*gamma^2) / 9 ) );
    
    % --------------------------------------------------------------------
    % solve the ADMM sub-problem for the primal variable u (linear system)
    % --------------------------------------------------------------------
        
        % store u_old for future computation of iterates relative change
        u_old = u; 
        % compute the right-hand side of the linear system
        rhs   = M .* (r - lambda_r / beta_r + b);
        % compute the new iterate
        u     = reshape( full( SUB_X_CM \ rhs(:) ) , h, w );
                   
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
        Su_b      = M .* u - b;
    
    % ------------------------------------------------------------------
    % update dual variables (Lagrange multipliers) lambdas (dual ascent)
    % ------------------------------------------------------------------
        
        % Lagrange multipliers associated with the auxiliary variable r
        lambda_r = lambda_r - beta_r * (r - Su_b);
        
    % --------------------------------------
    % (eventually) compute/store/show debugs
    % --------------------------------------
    if ( itrs_debugs == 1 )
        % compute and store iterations-based debug quantities
        Su_b_vec           = Su_b( M > 0 );
        J_regs(itrs+1)     = (1/2) * sum( Dhu(:).^2 + Dvu(:).^2 );
        J_fids(itrs+1)     = mu * sum( abs(Su_b_vec) );
        rel_chgs(itrs)     = rel_chg;
        Su_b_means(itrs+1) = mean( Su_b_vec );
        taus(itrs+1)       = sqrt( sum( Su_b_vec.^2 ) ) / (sqrt(d_M)*n_sigma);
        snrs(itrs+1)       = compute_snr(u,u_true);
        ssims(itrs+1)      = ssim(u,u_true);
        % Matlab command window debugs
        fprintf('\nIT %04d:  RC = %9.7f, J = %13.6f, tau* = %6.4f, [ISNR,ISSIM] = [%6.3f , %6.4f]',...
            itrs,rel_chgs(itrs),(J_regs(itrs+1)+J_fids(itrs+1)),taus(itrs+1),snrs(itrs+1)-snrs(1),ssims(itrs+1)-ssims(1));
    end
    
end

% one of the two stopping criteria satisfied --> end iterations  
% and, before returning the inpainted image u and the total number 
% of performed ADMM iterations itrs to the main script, show 
% (eventually if itrs_debugs = 1) some iterations-based debug quantities

% show iterations-based debug quantities
if ( itrs_debugs == 1 )

    % first, shrink arrays where we stored iterations-based debug quantities
    J_regs        = J_regs(1:itrs+1);
    J_fids        = J_fids(1:itrs+1);
    rel_chgs      = rel_chgs(1:itrs);
    Su_b_means    = Su_b_means(1:itrs+1);
    taus          = taus(1:itrs+1);
    snrs          = snrs(1:itrs+1);
    ssims         = ssims(1:itrs+1);

    % then show all debugs
    fprintf('\n');

    plots_x_min = 0;
    plots_x_max = itrs;

    figure('Name','Figure 30:  TIK-L1-U-ADMM:   ITERATIONS-BASED DEBUGS','NumberTitle','off');
    set(gcf,'Position',get(0,'ScreenSize'));

    subplot(2,3,1)
    plot(0:itrs,J_regs,'b');
    xlabel('k  (ADMM iteration index)');
    ylabel('$\mathrm{TIK}\left(u^{(k)}\right) = (1/2)\left\| D u \right\|_2^2$','interpreter','latex');
    plot1_y_min = min( J_regs(2:end) );
    plot1_y_max = max( J_regs(2:end) );
    plot1_y_min = plot1_y_min - 0.05 * (plot1_y_max - plot1_y_min);
    plot1_y_max = plot1_y_max + 0.05 * (plot1_y_max - plot1_y_min);
    axis([plots_x_min,plots_x_max,plot1_y_min,plot1_y_max]);

    subplot(2,3,2)
    plot(0:itrs,J_fids,'k');
    xlabel('k  (ADMM iteration index)');
    ylabel('$\mu \,\, \mathrm{L}_1\left(u^{(k)};b,S\right) = \mu \,\, \left\|Su^{(k)}-b\right\|_1$','interpreter','latex');
    plot2_y_min = min( J_fids(2:end) );
    plot2_y_max = max( J_fids(2:end) );
    plot2_y_min = plot2_y_min - 0.05 * (plot2_y_max - plot2_y_min);
    plot2_y_max = plot2_y_max + 0.05 * (plot2_y_max - plot2_y_min);
    axis([plots_x_min,plots_x_max,plot2_y_min,plot2_y_max]);

    subplot(2,3,3)
    plot(0:itrs,J_regs + J_fids,'r');
    xlabel('k  (ADMM iteration index)');
    ylabel('$J\left(u^{(k)};\mu\right) = \mathrm{TIK}\left(u^{(k)}\right) + \mu \, \mathrm{L}_1\left(u^{(k)};b,S\right)$','interpreter','latex');
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
    ylabel('$\tau^{(k)} = \left\|S u^{(k)} - b\right\|_2 / (\sqrt{d} \: \sigma_n)$','interpreter','latex');
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


