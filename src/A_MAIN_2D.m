% MAIN script of Matlab code for the course (see course slides):
%
%     MATHEMATICAL AND MACHINE LEARNING METHODS 
%               IN IMAGING (MOD.2) 
%
% By this code, one can run/test the 4 Unconstrained (U) variational models:
%
% TIK-L2-U, TIK-L1-U, TV-L2-U, TV-L1-U
%
% and the discrepancy-constrained model:
%
% TV-L2-DC 
%
% for reconstructing images corrupted by blurring/masking and noise 
% (that is, for solving the RESTORATION and INPAINTING / SUPER-RESOLUTION
% inverse problems), as selected and applied to the true uncorrupted image 
% in the files:
%
% B1_REST_GENERATE_DATA_2D.m (for RESTORATION)
% B2_INPT_GENERATE_DATA_2D.m (for INPAINTING / SUPER-RESOLUTION)
%
% In particular, as illustrated in the slides of the course, the 
% quadratic TIK-L2 model is solved in a direct efficient manner 
% by using the 2D Discrete Fourier Transform (DFT) implemented
% by the very efficient 2D Fast Fourier Transform (FFT) -> O(d log(d))
% for the restoration, by the (sparse) Cholesky method (in particular, 
% by using the backslash \ Matlab command) for inpainting / super-resolution.
%
% The other 3 unconstrained models, TIK-L1, TV-L2, TV-L1, are instead 
% solved by the state-of-the-art iterative minimization approach ADMM. 
%
% Finally, the Discrepancy-Constrained (DC) version of the TV-L2 model for 
% restoration can also be run/tested, solved numerically by the ADMM approach.
% 
% Author: 
% Alessandro Lanza, University of Bologna, Department of Mathematics, 
% alessandro.lanza2@unibo.it

% clear/initialize Matlab environment
close all; clear all; clc; format short;

% initialize Matlab pseudo-random numbers generators
rng('default'); rng(3);

% choose if "profiling" (that means "measuring") the processing 
% CPU times of the overall Matlab code we are going to run:
% set PROFILE_CPU_TIMES = 1 for performing the CPU time profiling
PROFILE_CPU_TIMES = 0;
if ( PROFILE_CPU_TIMES == 1 )
    profile off; % stop previously started profiling
    profile on;  % start a new profiling
end

% set the INVERSE PROBLEM to be tested:
% set 0 to test IMAGE RESTORATION, 1 to test IMAGE INPAINTING
  TEST_REST_OR_INPT = 1;

% set if running the fast version of ADMM algorithms for 
% INPAITING/SUPER-RESOLUTION (based on a suitable split)
  FAST_INPT = 0;

% set MODELS and ALGORITHMS to be tested, where we use notations
% U = 'UNCONSTRAINED', DC = 'DISCREPANCY-CONSTRAINED':
% set 1 to test, 0 not to test (you can set 1 for more than one approach).
    % Unconstrained quadratic TIK-L2 model by direct 2D FFT (REST) or PCG / DIRECT (INPT)
      TEST_TIK_L2_U_DIR   = 0;
    % Unconstrained TIK-L1, TV-L2, TV-L1, models by (iterative) ADMM approach
      TEST_TIK_L1_U_ADMM  = 0;
      TEST_TIK_L15_U_ADMM = 1;
      TEST_TV_L1_U_ADMM   = 0;
      TEST_TV_L15_U_ADMM  = 0;
      TEST_TV_L2_U_ADMM   = 0;
    % Discrepancy-Constrained TV-L2 model by (iterative) ADMM approach
      TEST_TV_L2_DC_ADMM  = 0;

% choose what to SHOW and/or to SAVE during the code execution
    % what to SHOW (in Matlab figures)
    SHOW_CORRUPTIONS                    = 1;
    SHOW_NOISE_HISTOGRAMS               = 0;
    SHOW_ALL_RESTORED_IMAGES            = 0;
    SHOW_BEST_RESTORED_IMAGES           = 1;
    SHOW_PERFORMANCE_GRAPHS_VS_MU       = 1;
    SHOW_TAU_STAR_VS_MU                 = 0; %not important
    SHOW_PERFORMANCE_GRAPHS_VS_TAU_STAR = 0; %not important
    % where (in which folder) and what (which images) to SAVE
    SAVE_FOLDER                         = './output';
    SAVE_ALL_RESTORED_IMAGES            = 0;
    SAVE_BEST_RESTORED_IMAGES           = 0; %suggested, only save the best restored images
    % scaling factor for visualization/saving of
    % the reconstruction absolute error images
    AEI_VIS_SF = 5; %don't change, not important, scaling factor to better visualize the absolute error images

% generate - and eventually visualize - the original 
% and corrupted (that is, blurred/masked and noisy) 
% image b by calling the GENERATE_DATA_2D files
if ( TEST_REST_OR_INPT == 0 )
    B1_REST_GENERATE_DATA_2D;
    Ax = @(x) real(ifft2( K_DFT .* fft2(x) ));
    INV_PROB_NAME_LONG  = 'RESTORATION';
    INV_PROB_NAME_SHORT = 'REST';
else
    B2_INPT_GENERATE_DATA_2D;
    Ax = @(x) M .* x; 
    INV_PROB_NAME_LONG  = 'INPAINTING';
    INV_PROB_NAME_SHORT = 'INPT';
end

% ------------------------------------------------------------------
% SET MODEL PARAMETERS:
% FOR ALL MODELS, CHOOSE WHICH VALUES OF THEIR FREE SCALAR PARAMETER
% TO TEST, NAMELY THE REGULARIZATION PARAMETER mu FOR UNCONSTRAINED 
% MODELS, THE DISCREPANCY PARAMETER tau FOR CONSTRAINED MODELS
% ------------------------------------------------------------------

% for UNCONSTRAINED MODELS: 
% choose the values of the regularization parameter mu to be tested
    % TIK_L2_U
      TIK_L2_U_mus_min  = 0.05;
      TIK_L2_U_mus_max  = 0.5;
      TIK_L2_U_mus_n    = 81;
      TIK_L2_U_mus      = linspace(TIK_L2_U_mus_min,TIK_L2_U_mus_max,TIK_L2_U_mus_n);

    % TIK_L1_U
      TIK_L1_U_mus_min  = 0.010;
      TIK_L1_U_mus_max  = 0.035;
      TIK_L1_U_mus_n    = 6;
      TIK_L1_U_mus      = linspace(TIK_L1_U_mus_min,TIK_L1_U_mus_max,TIK_L1_U_mus_n);

    % TIK_L15_U
      TIK_L15_U_mus_min  = 0.01;
      TIK_L15_U_mus_max  = 10;
      TIK_L15_U_mus_n    = 25;
      TIK_L15_U_mus      = linspace(TIK_L15_U_mus_min,TIK_L15_U_mus_max,TIK_L15_U_mus_n);

    % TV_L1_U
      TV_L1_U_mus_min   = 0.5;  %0.5;
      TV_L1_U_mus_max   = 4;    %4.0;
      TV_L1_U_mus_n     = 2;
      TV_L1_U_mus       = linspace(TV_L1_U_mus_min,TV_L1_U_mus_max,TV_L1_U_mus_n); 

    % TV_L15_U
      TV_L15_U_mus_min   = 0.05;  %0.5;
      TV_L15_U_mus_max   = 30;    %4.0;
      TV_L15_U_mus_n     = 15;
      TV_L15_U_mus       = linspace(TV_L15_U_mus_min,TV_L15_U_mus_max,TV_L15_U_mus_n);

    % TV_L2_U
      TV_L2_U_mus_min   = 80;
      TV_L2_U_mus_max   = 200;
      TV_L2_U_mus_n     = 6;
      TV_L2_U_mus       = linspace(TV_L2_U_mus_min,TV_L2_U_mus_max,TV_L2_U_mus_n);    
  
% for CONSTRAINED MODELS: 
% choose the values of the discrepancy parameter tau to be tested
    % TV_L2_DC
      TV_L2_DC_taus_min = 1;
      TV_L2_DC_taus_max = 1;
      TV_L2_DC_taus_n   = 1;
      TV_L2_DC_taus     = linspace(TV_L2_DC_taus_min,TV_L2_DC_taus_max,TV_L2_DC_taus_n);

% ------------------------------------------------------------------
% SET MINIMIZATION ALGORITHMS PARAMETERS:
% FOR ALL MINIMIZATION ALGORITHMS (ADMM and QMM) APPLIED TO ALL 
% MODELS, CHOOSE THE VALUES OF ALL THE PARAMETERS WHICH DEFINE 
% COMPLETELY HOW THE ALGORITHMS ARE APPLIED FOR THE SOLUTION
% ------------------------------------------------------------------
    
% TIK-L2-U-DIR:
% no algorithm parameter to be set, as the solution is
% computed directly (no iterations) by using the 2D FFT
% for RESTORATION, \ (backslash) for INPAINTING

% common parameters for all iterative algorithms (ADMM and ...):
% initial iterate and iterations stopping criteria
  u0              = b;       % initial iterate
  itrs_th         = 1000;    % stopping criterion 1: threshold (maximum number) for iterations
  itrs_rel_chg_th = 10^(-4.5); % stopping criterion 2: threshold for iterates relative change 
        
% parameters of the ADMM iterative approach applied to solving the 
% unconstrained TIK-L1-U, TIK-L15-U, TV-L1-U, TV-L15-U and TV-L2-U models 
% and the discrepancy-constrained TV-L2-DC model (see course slides)
    % TIK-L1-U
      TIK_L1_U_ADMM_beta_r = 0.2;

    % TIK-L15-U
      TIK_L15_U_ADMM_beta_r = 0.2;

    % TV-L1-U
      TV_L1_U_ADMM_beta_t  = 3;
      TV_L1_U_ADMM_beta_r  = 3;

    % TV-L15-U
      TV_L15_U_ADMM_beta_t  = 3;
      TV_L15_U_ADMM_beta_r  = 2;
    
    % TV-L2-U
      TV_L2_U_ADMM_beta_t  = 1;

    % TV-L2-DC
      TV_L2_DC_ADMM_beta_t = 50;
      TV_L2_DC_ADMM_beta_r = 50;
    
% ------------------------------------------------------------------
% ALLOCATE (CREATE AND INITIALIZE TO ZERO) THE DATA STRUCTURES (ARRAYS)
% WHERE WE WILL STORE SOME IMPORTANT SCALAR MEASURES OF PERFORMANCE 
% FOR ALL THE SELECTED MODELS, ALL THE SELECTED VALUES OF THEIR FREE 
% SCALAR PARAMETER (mu or tau), AND ALL THE SELECTED MINIMIZATION 
% ALGORITHMS (DIRECT or ITERATIVE ADMM or ...).
% THE SCALAR MEASURES OF PERFORMANCE ARE:
% ACC.1)  ISNR  OF THE RESTORED IMAGE (MEASURE OF RESTORATION ACCURACY);
% ACC.2)  ISSIM OF THE RESTORED IMAGE (MEASURE OF RESTORATION ACCURACY);
% EFF.1)  NUMBER OF PERFORMED ALGORITHM ITERATIONS (MEASURE OF RESTORATION EFFICIENCY)
% EFF.2)  CPU TIME SPENT BY THE ALGORITHMS (MEASURE OF RESTORATION EFFICIENCY)
% INSP.1) DISCREPANCY ASSOCIATED TO THE RESTORED IMAGE RESIDUAL (INSPECTION)
% ------------------------------------------------------------------

if (TEST_TIK_L2_U_DIR == 1)
    TIK_L2_U_DIR_ISNRS      = zeros(1,TIK_L2_U_mus_n);
    TIK_L2_U_DIR_ISSIMS     = zeros(1,TIK_L2_U_mus_n);
    TIK_L2_U_DIR_ITRS       =  ones(1,TIK_L2_U_mus_n);
    TIK_L2_U_DIR_CPU_TIMES  = zeros(1,TV_L2_U_mus_n);
    TIK_L2_U_DIR_TAU_STARS  = zeros(1,TIK_L2_U_mus_n);
end
%
if (TEST_TIK_L1_U_ADMM == 1)
    TIK_L1_U_ADMM_ISNRS     = zeros(1,TIK_L1_U_mus_n);
    TIK_L1_U_ADMM_ISSIMS    = zeros(1,TIK_L1_U_mus_n);
    TIK_L1_U_ADMM_ITRS      = zeros(1,TIK_L1_U_mus_n);
    TIK_L1_U_ADMM_CPU_TIMES = zeros(1,TIK_L1_U_mus_n);
    TIK_L1_U_ADMM_TAU_STARS = zeros(1,TIK_L1_U_mus_n);
end
%
if (TEST_TIK_L15_U_ADMM == 1)
    TIK_L15_U_ADMM_ISNRS     = zeros(1,TIK_L15_U_mus_n);
    TIK_L15_U_ADMM_ISSIMS    = zeros(1,TIK_L15_U_mus_n);
    TIK_L15_U_ADMM_ITRS      = zeros(1,TIK_L15_U_mus_n);
    TIK_L15_U_ADMM_CPU_TIMES = zeros(1,TIK_L15_U_mus_n);
    TIK_L15_U_ADMM_TAU_STARS = zeros(1,TIK_L15_U_mus_n);
end
%
if (TEST_TV_L1_U_ADMM == 1)
    TV_L1_U_ADMM_ISNRS      = zeros(1,TV_L1_U_mus_n);
    TV_L1_U_ADMM_ISSIMS     = zeros(1,TV_L1_U_mus_n);
    TV_L1_U_ADMM_ITRS       = zeros(1,TV_L1_U_mus_n);
    TV_L1_U_ADMM_CPU_TIMES  = zeros(1,TV_L1_U_mus_n);
    TV_L1_U_ADMM_TAU_STARS  = zeros(1,TV_L1_U_mus_n);
end
if (TEST_TV_L15_U_ADMM == 1)
    TV_L15_U_ADMM_ISNRS      = zeros(1,TV_L15_U_mus_n);
    TV_L15_U_ADMM_ISSIMS     = zeros(1,TV_L15_U_mus_n);
    TV_L15_U_ADMM_ITRS       = zeros(1,TV_L15_U_mus_n);
    TV_L15_U_ADMM_CPU_TIMES  = zeros(1,TV_L15_U_mus_n);
    TV_L15_U_ADMM_TAU_STARS  = zeros(1,TV_L15_U_mus_n);
end
%
if (TEST_TV_L2_U_ADMM == 1)
    TV_L2_U_ADMM_ISNRS      = zeros(1,TV_L2_U_mus_n);
    TV_L2_U_ADMM_ISSIMS     = zeros(1,TV_L2_U_mus_n);
    TV_L2_U_ADMM_ITRS       = zeros(1,TV_L2_U_mus_n);
    TV_L2_U_ADMM_CPU_TIMES  = zeros(1,TV_L2_U_mus_n);
    TV_L2_U_ADMM_TAU_STARS  = zeros(1,TV_L2_U_mus_n);
end
%
if (TEST_TV_L2_DC_ADMM == 1)
    TV_L2_DC_ADMM_ISNRS     = zeros(1,TV_L2_DC_taus_n);
    TV_L2_DC_ADMM_ISSIMS    = zeros(1,TV_L2_DC_taus_n);
    TV_L2_DC_ADMM_ITRS      = zeros(1,TV_L2_DC_taus_n);
    TV_L2_DC_ADMM_CPU_TIMES = zeros(1,TV_L2_DC_taus_n);
    TV_L2_DC_ADMM_TAU_STARS = zeros(1,TV_L2_DC_taus_n);
end

  
% in order to compute and store all the ISNR and ISSIM values above,
% we will need the SNR and the SSIM of the initial iterate u0, that 
% are the same for all models/algorithms, hence we compute them here
% once and for all. 
  u0_SNR      = compute_snr(u0,u_true);
  u0_SSIM     = ssim(u0,u_true);


% ------------------------------------------------------------------
% NOW WE CAN RUN THE CODE FOR TESTING ALL THE SELECTED VARIATIONAL 
% MODELS, WITH ALL THE SELECTED VALUES OF THEIR FREE SCALAR PARAMETER, 
% SOLVED NUMERICALLY BY ALL THE SELECTED MINIMIZATION ALGORITHMS
% ------------------------------------------------------------------

% first of all, independently of selected models/algorithms, 
% if one of the two SAVE_... is activated, save the true, the 
% blurred/masked and the observed (blurred/masked + noisy) images 
% if ( ( SAVE_ALL_RESTORED_IMAGES == 1 ) || ( SAVE_BEST_RESTORED_IMAGES == 1 ) )
%     imwrite(uint8(IMS_VIS_SF*u_true),sprintf('%s/%s_IM%02d_A_ORIG.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id));
%     if ( TEST_REST_OR_INPT == 0 )
%         imwrite(uint8(IMS_VIS_SF*Au_true),sprintf('%s/%s_IM%02d_B_Au.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id));
%         imwrite(uint8(IMS_VIS_SF*b),sprintf('%s/%s_IM%02d_C_CORR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id));
%         imwrite(uint8( 127.5+IMS_VIS_SF*(b - Au_true) ),sprintf('%s/%s_IM%02d_C_NOISE.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id));
%     else
%         imwrite(uint8( MASK_IMAGE_COL(IMS_VIS_SF*Au_true,M,M_COL) ),sprintf('%s/%s_IM%02d_B_Au.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id));
%         imwrite(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ),sprintf('%s/%s_IM%02d_C_CORR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id));
%         imwrite(uint8( 255 * M ),sprintf('%s/%s_IM%02d_C_MASK.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id));
%         imwrite(uint8( MASK_IMAGE_COL(127.5+IMS_VIS_SF*(b - Au_true),M,M_COL) ),sprintf('%s/%s_IM%02d_C_NOISE_MASKED.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id));
%         imwrite(uint8( 127.5+IMS_VIS_SF*(b - Au_true) ),sprintf('%s/%s_IM%02d_C_NOISE.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id));
%     end
% end

%--------------------------------------------------------------------------
   
% TEST TIK_L2_U MODEL SOLVED DIRECTLY BY 2D DFT (REST) OR CHOLESKY (INPT)
if ( TEST_TIK_L2_U_DIR == 1 )
    % compute (in order to print it) tau_star of u0
    Au0_b       = Ax(u0) - b;
    u0_tau_star = norm( Au0_b(:) ) / (sqrt(d)*n_sigma);
    % print preliminaries of model/algorithm/initial iterate in the command window
    fprintf('\n\n\nIMAGE %s by TIK-L2-U MODEL solved directly:\n',INV_PROB_NAME_LONG);
    fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);
       
    % cycle over all the selected different mu values
    TIK_L2_U_DIR_BEST_ISNR  = -10^5; 
    TIK_L2_U_DIR_BEST_ISSIM = -10^5;
    for mu_i = 1:TIK_L2_U_mus_n
        mu = TIK_L2_U_mus(mu_i);
        % for the current value of the regularization parameter mu, 
        % compute the reconstructed image u_star(mu) by direct (that is, 
        % non-iterative) solution of the 1st-order optimality conditions 
        % (linear system) of the quadratic unconstrained TIK-L2 model 
        % using 2D DFT (restoration) or Cholesky (inpainting)
        t0 = cputime;
        
        if ( TEST_REST_OR_INPT == 0 )
            u_star = REST_TIK_L2_U_DIR(b,blur_k,mu);
        else
            u_star = INPT_TIK_L2_U_DIR(b,M,mu);
        end
        
        % compute/store the current reconstruction scalar performance measures
        TIK_L2_U_DIR_CPU_TIMES(mu_i) = cputime - t0;
        TIK_L2_U_DIR_ISNRS(mu_i)     = compute_snr(u_star,u_true) - u0_SNR;
        TIK_L2_U_DIR_ISSIMS(mu_i)    = ssim(u_star,u_true) - u0_SSIM;
        Au_star_b                    = Ax(u_star) - b;
        TIK_L2_U_DIR_TAU_STARS(mu_i) = norm( Au_star_b(:) ) / (sqrt(d)*n_sigma);
        % update the best reconstruction results in terms of accuracy measures ISNR and SSIM 
        if ( TIK_L2_U_DIR_ISNRS(mu_i) > TIK_L2_U_DIR_BEST_ISNR )
            TIK_L2_U_DIR_BEST_ISNR         = TIK_L2_U_DIR_ISNRS(mu_i);
            TIK_L2_U_DIR_BEST_ISNR_u_star  = u_star;
            TIK_L2_U_DIR_BEST_ISNR_mu_i    = mu_i;
        end
        if ( TIK_L2_U_DIR_ISSIMS(mu_i) > TIK_L2_U_DIR_BEST_ISSIM )
            TIK_L2_U_DIR_BEST_ISSIM        = TIK_L2_U_DIR_ISSIMS(mu_i);
            TIK_L2_U_DIR_BEST_ISSIM_u_star = u_star;
            TIK_L2_U_DIR_BEST_ISSIM_mu_i   = mu_i;
        end
        % print in the Matlab command window the current reconstruction results
        fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,  [ISNR,ISSIM] = [%6.3f , %6.4f]',...
            mu_i,TIK_L2_U_mus_n,mu,1,TIK_L2_U_DIR_TAU_STARS(mu_i),TIK_L2_U_DIR_ISNRS(mu_i),TIK_L2_U_DIR_ISSIMS(mu_i));
        % show in a Matlab figure the current reconstruction results
        if ( SHOW_ALL_RESTORED_IMAGES == 1 )
            if ( mu_i == 1 )
                figure('Name','Figure 11:  TIK-L2-U-DIR:   mu-DEPENDENT RESULTS','NumberTitle','off');
                set(gcf,'Position',get(0,'ScreenSize'));
            end
            subplot(2,2,1)
            imshow(uint8(IMS_VIS_SF*u_true));
            title(sprintf('IM%02d (%d x %d):   u_{true}',im_id,h,w));
            subplot(2,2,3)
            if ( TEST_REST_OR_INPT == 0 )
                imshow(uint8(IMS_VIS_SF*b));
            else
                imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
            end
            title(sprintf('b = N(A u_{true}):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM));
            subplot(2,2,2)
            imshow(uint8(IMS_VIS_SF*u_star));
            title(sprintf('u^*(mu = %7.2f):  [ISNR,ISSIM] = [%6.3f,%6.4f]',mu,TIK_L2_U_DIR_ISNRS(mu_i),TIK_L2_U_DIR_ISSIMS(mu_i)));
            subplot(2,2,4)
            imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star) )); 
            title(sprintf('| u^*(mu = %7.2f) - u_{true} |',mu));
            pause(0.1);
        end
        % save on the hard disk the current reconstruction results
        if ( SAVE_ALL_RESTORED_IMAGES == 1 )
            imwrite(uint8(IMS_VIS_SF*u_star),sprintf('%s/%s_IM%02d_D1_TIK_L2_U_DIR_A_REST_mu_%07.2f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,mu));
            imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star)),sprintf('%s/%s_IM%02d_D1_TIK_L2_U_DIR_B_AERR_mu_%07.2f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,mu));
        end
    end
    % print in the Matlab command window the best reconstruction results in terms of ISNR and ISSIM
    TIK_L2_U_DIR_BEST_ISNR_mu        = TIK_L2_U_mus(TIK_L2_U_DIR_BEST_ISNR_mu_i);
    TIK_L2_U_DIR_BEST_ISSIM_mu       = TIK_L2_U_mus(TIK_L2_U_DIR_BEST_ISSIM_mu_i);
    TIK_L2_U_DIR_BEST_ISNR_itrs      = 1;%TIK_L2_U_DIR_ITRS(TIK_L2_U_DIR_BEST_ISNR_mu_i);
    TIK_L2_U_DIR_BEST_ISSIM_itrs     = 1;%TIK_L2_U_DIR_ITRS(TIK_L2_U_DIR_BEST_ISSIM_mu_i);
    TIK_L2_U_DIR_BEST_ISNR_tau_star  = TIK_L2_U_DIR_TAU_STARS(TIK_L2_U_DIR_BEST_ISNR_mu_i);
    TIK_L2_U_DIR_BEST_ISSIM_tau_star = TIK_L2_U_DIR_TAU_STARS(TIK_L2_U_DIR_BEST_ISSIM_mu_i);        
    fprintf('\nBEST ISNR and ISSIM RESULTS:');
    fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,   BEST ISNR   =  %6.3f',...
        TIK_L2_U_DIR_BEST_ISNR_mu_i,TIK_L2_U_mus_n,TIK_L2_U_DIR_BEST_ISNR_mu,...
        TIK_L2_U_DIR_BEST_ISNR_itrs,TIK_L2_U_DIR_BEST_ISNR_tau_star,TIK_L2_U_DIR_BEST_ISNR);
    fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,   BEST ISSIM  =   %6.4f',...
        TIK_L2_U_DIR_BEST_ISSIM_mu_i,TIK_L2_U_mus_n,TIK_L2_U_DIR_BEST_ISSIM_mu,...
        TIK_L2_U_DIR_BEST_ISSIM_itrs,TIK_L2_U_DIR_BEST_ISSIM_tau_star,TIK_L2_U_DIR_BEST_ISSIM);
    % show in a Matlab figure the best reconstruction results in terms of ISNR and ISSIM
    if ( SHOW_BEST_RESTORED_IMAGES == 1 )
        figure('Name','Figure 12:  TIK-L2-U-DIR:  BEST RESULTS in terms of ISNR and ISSIM','NumberTitle','off');
        set(gcf,'Position',get(0,'ScreenSize'));
        subplot(2,3,1)
        imshow(uint8(IMS_VIS_SF*u_true));
        title(sprintf('IM%02d (%d x %d):   u_{true}',im_id,h,w));
        subplot(2,3,4)
        if ( TEST_REST_OR_INPT == 0 )
            imshow(uint8(IMS_VIS_SF*b));
        else
            imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
        end
        title(sprintf('b = N(A u_{true}):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM));
        subplot(2,3,2)
        imshow(uint8(IMS_VIS_SF*TIK_L2_U_DIR_BEST_ISNR_u_star));
        title(sprintf('u^*(mu = %7.2f):  BEST ISNR = %6.3f',TIK_L2_U_DIR_BEST_ISNR_mu,TIK_L2_U_DIR_BEST_ISNR));
        subplot(2,3,5)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TIK_L2_U_DIR_BEST_ISNR_u_star) )); 
        title(sprintf('| u^*(mu = %7.2f) - u_{true} |',TIK_L2_U_DIR_BEST_ISNR_mu));
        subplot(2,3,3)
        imshow(uint8(IMS_VIS_SF*TIK_L2_U_DIR_BEST_ISSIM_u_star));
        title(sprintf('u^*(mu = %7.2f):  BEST ISSIM = %6.4f',TIK_L2_U_DIR_BEST_ISSIM_mu,TIK_L2_U_DIR_BEST_ISSIM));
        subplot(2,3,6)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TIK_L2_U_DIR_BEST_ISSIM_u_star) )); 
        title(sprintf('| u^*(mu = %7.2f) - u_{true} |',TIK_L2_U_DIR_BEST_ISSIM_mu));
    end
    % save on the hard disk the best reconstruction results in terms of ISNR and ISSIM
    if ( SAVE_BEST_RESTORED_IMAGES == 1 )
        imwrite(uint8(IMS_VIS_SF*TIK_L2_U_DIR_BEST_ISNR_u_star),sprintf('%s/%s_IM%02d_D1_TIK_L2_U_DIR_A_REST_mu_%07.2f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TIK_L2_U_DIR_BEST_ISNR_mu));
        imwrite(uint8(IMS_VIS_SF*TIK_L2_U_DIR_BEST_ISSIM_u_star),sprintf('%s/%s_IM%02d_D1_TIK_L2_U_DIR_A_REST_mu_%07.2f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TIK_L2_U_DIR_BEST_ISSIM_mu));
        imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TIK_L2_U_DIR_BEST_ISNR_u_star)),sprintf('%s/%s_IM%02d_D1_TIK_L2_U_DIR_B_AERR_mu_%07.2f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TIK_L2_U_DIR_BEST_ISNR_mu));
        imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TIK_L2_U_DIR_BEST_ISSIM_u_star)),sprintf('%s/%s_IM%02d_D1_TIK_L2_U_DIR_B_AERR_mu_%07.2f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TIK_L2_U_DIR_BEST_ISSIM_mu));
    end
        
    % show in Matlab figures the performance graphs
    % show...only if we have selected at least two different mu values
    mus_min = min(TIK_L2_U_mus); mus_max = max(TIK_L2_U_mus);
    if ( mus_max > mus_min )
        if ( SHOW_PERFORMANCE_GRAPHS_VS_MU == 1 )
            figure('Name','Figure 13:  TIK-L2-U-DIR:  PERFORMANCE GRAPHS  versus  MU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TIK-L2-U-DIR:    ACCURACY MEASURES  versus  \mu');
            yyaxis left;
            plot_y_min = min(TIK_L2_U_DIR_ISNRS); plot_y_max = max(TIK_L2_U_DIR_ISNRS); 
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L2_U_mus,TIK_L2_U_DIR_ISNRS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ISNR(\mu)  [dB]'); 
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TIK_L2_U_DIR_ISSIMS); plot_y_max = max(TIK_L2_U_DIR_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L2_U_mus,TIK_L2_U_DIR_ISSIMS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ISSIM(\mu)');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TIK-L2-U-DIR:    EFFICIENCY MEASURES  versus  \mu');
            yyaxis left;
            plot_y_min = min(TIK_L2_U_DIR_ITRS); plot_y_max = max(TIK_L2_U_DIR_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L2_U_mus,TIK_L2_U_DIR_ITRS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ITERS(\mu)');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TIK_L2_U_DIR_CPU_TIMES); plot_y_max = max(TIK_L2_U_DIR_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L2_U_mus,TIK_L2_U_DIR_CPU_TIMES,'-o','markersize',3);
            xlabel('\mu'); ylabel('CPU-TIME(\mu)  [secs]');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_TAU_STAR_VS_MU == 1 )
            figure('Name','Figure 14:  TIK-L2-U-DIR:  TAU*  versus  MU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));            
            plot_y_min = min(TIK_L2_U_DIR_TAU_STARS); plot_y_max = max(TIK_L2_U_DIR_TAU_STARS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L2_U_mus,TIK_L2_U_DIR_TAU_STARS,'-om','markersize',3);
            xlabel('\mu'); ylabel('\tau^*(\mu)');
            title('TIK-L2-U-DIR:    TAU*  versus  \mu');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_PERFORMANCE_GRAPHS_VS_TAU_STAR == 1 )
            tau_stars_min = min(TIK_L2_U_DIR_TAU_STARS); tau_stars_max = max(TIK_L2_U_DIR_TAU_STARS);
            figure('Name','Figure 15:  TIK-L2-U-DIR:  PERFORMANCE GRAPHS  versus  TAU*(MU)','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TIK-L2-U-DIR:    ACCURACY MEASURES  versus  \tau^*(\mu)');
            yyaxis left;
            plot_y_min = min(TIK_L2_U_DIR_ISNRS); plot_y_max = max(TIK_L2_U_DIR_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L2_U_DIR_TAU_STARS,TIK_L2_U_DIR_ISNRS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ISNR(\tau^*(\mu))  [dB]'); 
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TIK_L2_U_DIR_ISSIMS); plot_y_max = max(TIK_L2_U_DIR_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L2_U_DIR_TAU_STARS,TIK_L2_U_DIR_ISSIMS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ISSIM(\tau^*(\mu))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TIK-L2-U-DIR:    EFFICIENCY MEASURES  versus  \tau^*(\mu)');
            yyaxis left;
            plot_y_min = min(TIK_L2_U_DIR_ITRS); plot_y_max = max(TIK_L2_U_DIR_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L2_U_DIR_TAU_STARS,TIK_L2_U_DIR_ITRS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ITERS(\tau^*(\mu))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TIK_L2_U_DIR_CPU_TIMES); plot_y_max = max(TIK_L2_U_DIR_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L2_U_DIR_TAU_STARS,TIK_L2_U_DIR_CPU_TIMES,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('CPU-TIME(\tau^*(\mu))  [secs]');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
        end
    end
end

%--------------------------------------------------------------------------

if ( TEST_TIK_L1_U_ADMM == 1 )
    % compute (in order to print it) tau_star of u0
    Au0_b       = Ax(u0) - b;
    u0_tau_star = norm( Au0_b(:) ) / (sqrt(d)*n_sigma);
    % print preliminaries of model/algorithm/initial iterate in the command window
    fprintf('\n\n\nIMAGE %s by TIK-L1-U MODEL solved by ADMM:\n',INV_PROB_NAME_LONG);
    if ( TIK_L1_U_mus_n == 1 ), itrs_debugs = 1; else, itrs_debugs = 0; 
        fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);    
    end
    % cycle over all the selected different mu values
    TIK_L1_U_ADMM_BEST_ISNR  = -10^5; 
    TIK_L1_U_ADMM_BEST_ISSIM = -10^5;
    for mu_i = 1:TIK_L1_U_mus_n
        mu = TIK_L1_U_mus(mu_i);
        % for the current value of the regularization parameter mu,  
        % compute the reconstructed image u_star(mu) by solving the model by ADMM        
        t0 = cputime;
        
        if ( TEST_REST_OR_INPT == 0 )
            [u_star,itrs] = REST_TIK_L1_U_ADMM(b,blur_k,mu,TIK_L1_U_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
        else
            if ( FAST_INPT == 1 )
                [u_star,itrs] = INPT_TIK_L1_U_ADMM_FAST(b,M,mu,TIK_L1_U_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
            else
                [u_star,itrs] = INPT_TIK_L1_U_ADMM(b,M,mu,TIK_L1_U_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
            end
        end
        % compute/store the current reconstruction scalar performance measures
        TIK_L1_U_ADMM_CPU_TIMES(mu_i) = cputime - t0;
        TIK_L1_U_ADMM_ISNRS(mu_i)     = compute_snr(u_star,u_true) - u0_SNR;
        TIK_L1_U_ADMM_ISSIMS(mu_i)    = ssim(u_star,u_true) - u0_SSIM;
        TIK_L1_U_ADMM_ITRS(mu_i)      = itrs;
        Au_star_b                     = Ax(u_star) - b;
        TIK_L1_U_ADMM_TAU_STARS(mu_i) = norm( Au_star_b(:) ) / (sqrt(d)*n_sigma);
        % update the best reconstruction results in terms of accuracy measures ISNR and SSIM 
        if ( TIK_L1_U_ADMM_ISNRS(mu_i) > TIK_L1_U_ADMM_BEST_ISNR )
            TIK_L1_U_ADMM_BEST_ISNR         = TIK_L1_U_ADMM_ISNRS(mu_i);
            TIK_L1_U_ADMM_BEST_ISNR_u_star  = u_star;
            TIK_L1_U_ADMM_BEST_ISNR_mu_i    = mu_i;
        end
        if ( TIK_L1_U_ADMM_ISSIMS(mu_i) > TIK_L1_U_ADMM_BEST_ISSIM )
            TIK_L1_U_ADMM_BEST_ISSIM        = TIK_L1_U_ADMM_ISSIMS(mu_i);
            TIK_L1_U_ADMM_BEST_ISSIM_u_star = u_star;
            TIK_L1_U_ADMM_BEST_ISSIM_mu_i   = mu_i;
        end
        % print in the Matlab command window the current reconstruction results
        if ( TIK_L1_U_mus_n == 1 )
            fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);
        end
        fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,  [ISNR,ISSIM] = [%6.3f , %6.4f]',...
            mu_i,TIK_L1_U_mus_n,mu,itrs,TIK_L1_U_ADMM_TAU_STARS(mu_i),TIK_L1_U_ADMM_ISNRS(mu_i),TIK_L1_U_ADMM_ISSIMS(mu_i));
        % show in a Matlab figure the current reconstruction results
        if ( SHOW_ALL_RESTORED_IMAGES == 1 )
            if ( mu_i == 1 )
                figure('Name','Figure 31:  TIK-L1-U-ADMM:   mu-DEPENDENT RESULTS','NumberTitle','off');
                set(gcf,'Position',get(0,'ScreenSize'));
            end
            subplot(2,2,1)
            imshow(uint8(IMS_VIS_SF*u_true));
            title(sprintf('IM%02d (%d x %d):   u_{true}',im_id,h,w));
            subplot(2,2,3)
            if ( TEST_REST_OR_INPT == 0 )
                imshow(uint8(IMS_VIS_SF*b));
            else
                imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
            end
            title(sprintf('b = N(A u_{true}):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM));
            subplot(2,2,2)
            imshow(uint8(IMS_VIS_SF*u_star));
            title(sprintf('u^*(mu = %7.2f):  [ISNR,ISSIM] = [%6.3f,%6.4f]',mu,TIK_L1_U_ADMM_ISNRS(mu_i),TIK_L1_U_ADMM_ISSIMS(mu_i)));
            subplot(2,2,4)
            imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star) )); 
            title(sprintf('| u^*(mu = %7.2f) - u_{true} |',mu));
            pause(0.1);
        end
        % save on the hard disk the current reconstruction results
        if ( SAVE_ALL_RESTORED_IMAGES == 1 )
            imwrite(uint8(IMS_VIS_SF*u_star),sprintf('%s/%s_IM%02d_D3_TIK_L1_U_ADMM_A_REST_mu_%07.2f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,mu));
            imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star)),sprintf('%s/%s_IM%02d_D3_TIK_L1_U_ADMM_B_AERR_mu_%07.2f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,mu));
        end
    end
    % print in the Matlab command window the best reconstruction results in terms of ISNR and ISSIM
    TIK_L1_U_ADMM_BEST_ISNR_mu        = TIK_L1_U_mus(TIK_L1_U_ADMM_BEST_ISNR_mu_i);
    TIK_L1_U_ADMM_BEST_ISSIM_mu       = TIK_L1_U_mus(TIK_L1_U_ADMM_BEST_ISSIM_mu_i);
    TIK_L1_U_ADMM_BEST_ISNR_itrs      = TIK_L1_U_ADMM_ITRS(TIK_L1_U_ADMM_BEST_ISNR_mu_i);
    TIK_L1_U_ADMM_BEST_ISSIM_itrs     = TIK_L1_U_ADMM_ITRS(TIK_L1_U_ADMM_BEST_ISSIM_mu_i);
    TIK_L1_U_ADMM_BEST_ISNR_tau_star  = TIK_L1_U_ADMM_TAU_STARS(TIK_L1_U_ADMM_BEST_ISNR_mu_i);
    TIK_L1_U_ADMM_BEST_ISSIM_tau_star = TIK_L1_U_ADMM_TAU_STARS(TIK_L1_U_ADMM_BEST_ISSIM_mu_i);        
    fprintf('\nBEST ISNR and ISSIM RESULTS:');
    fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,   BEST ISNR   =  %6.3f',...
        TIK_L1_U_ADMM_BEST_ISNR_mu_i,TIK_L1_U_mus_n,TIK_L1_U_ADMM_BEST_ISNR_mu,...
        TIK_L1_U_ADMM_BEST_ISNR_itrs,TIK_L1_U_ADMM_BEST_ISNR_tau_star,TIK_L1_U_ADMM_BEST_ISNR);
    fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,   BEST ISSIM  =   %6.4f',...
        TIK_L1_U_ADMM_BEST_ISSIM_mu_i,TIK_L1_U_mus_n,TIK_L1_U_ADMM_BEST_ISSIM_mu,...
        TIK_L1_U_ADMM_BEST_ISSIM_itrs,TIK_L1_U_ADMM_BEST_ISSIM_tau_star,TIK_L1_U_ADMM_BEST_ISSIM);
    % show in a Matlab figure the best reconstruction results in terms of ISNR and ISSIM
    if ( SHOW_BEST_RESTORED_IMAGES == 1 )
        figure('Name','Figure 32:  TIK-L1-U-ADMM:  BEST RESULTS in terms of ISNR and ISSIM','NumberTitle','off');
        set(gcf,'Position',get(0,'ScreenSize'));
        subplot(2,3,1)
        imshow(uint8(IMS_VIS_SF*u_true));
        title(sprintf('IM%02d (%d x %d):   u_{true}',im_id,h,w));
        subplot(2,3,4)
        if ( TEST_REST_OR_INPT == 0 )
            imshow(uint8(IMS_VIS_SF*b));
        else
            imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
        end
        title(sprintf('b = N(A u_{true}):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM));
        subplot(2,3,2)
        imshow(uint8(IMS_VIS_SF*TIK_L1_U_ADMM_BEST_ISNR_u_star));
        title(sprintf('u^*(mu = %7.2f):  BEST ISNR = %6.3f',TIK_L1_U_ADMM_BEST_ISNR_mu,TIK_L1_U_ADMM_BEST_ISNR));
        subplot(2,3,5)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TIK_L1_U_ADMM_BEST_ISNR_u_star) )); 
        title(sprintf('| u^*(mu = %7.2f) - u_{true} |',TIK_L1_U_ADMM_BEST_ISNR_mu));
        subplot(2,3,3)
        imshow(uint8(IMS_VIS_SF*TIK_L1_U_ADMM_BEST_ISSIM_u_star));
        title(sprintf('u^*(mu = %7.2f):  BEST ISSIM = %6.4f',TIK_L1_U_ADMM_BEST_ISSIM_mu,TIK_L1_U_ADMM_BEST_ISSIM));
        subplot(2,3,6)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TIK_L1_U_ADMM_BEST_ISSIM_u_star) )); 
        title(sprintf('| u^*(mu = %7.2f) - u_{true} |',TIK_L1_U_ADMM_BEST_ISSIM_mu));
    end
    % save on the hard disk the best reconstruction results in terms of ISNR and ISSIM
    if ( SAVE_BEST_RESTORED_IMAGES == 1 )
        imwrite(uint8(IMS_VIS_SF*TIK_L1_U_ADMM_BEST_ISNR_u_star),sprintf('%s/%s_IM%02d_D3_TIK_L1_U_ADMM_A_REST_mu_%07.2f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TIK_L1_U_ADMM_BEST_ISNR_mu));
        imwrite(uint8(IMS_VIS_SF*TIK_L1_U_ADMM_BEST_ISSIM_u_star),sprintf('%s/%s_IM%02d_D3_TIK_L1_U_ADMM_A_REST_mu_%07.2f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TIK_L1_U_ADMM_BEST_ISSIM_mu));
        imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TIK_L1_U_ADMM_BEST_ISNR_u_star)),sprintf('%s/%s_IM%02d_D3_TIK_L1_U_ADMM_B_AERR_mu_%07.2f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TIK_L1_U_ADMM_BEST_ISNR_mu));
        imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TIK_L1_U_ADMM_BEST_ISSIM_u_star)),sprintf('%s/%s_IM%02d_D3_TIK_L1_U_ADMM_B_AERR_mu_%07.2f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TIK_L1_U_ADMM_BEST_ISSIM_mu));
    end
    % show in a Matlab figure the performance graphs
    % show...only if we have selected at least two different mu values
    mus_min = min(TIK_L1_U_mus); mus_max = max(TIK_L1_U_mus);
    if ( mus_max > mus_min )
        if ( SHOW_PERFORMANCE_GRAPHS_VS_MU == 1 )
            figure('Name','Figure 33:  TIK-L1-U-ADMM:  PERFORMANCE GRAPHS  versus  MU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TIK-L1-U-ADMM:    ACCURACY MEASURES  versus  \mu');
            yyaxis left;
            plot_y_min = min(TIK_L1_U_ADMM_ISNRS); plot_y_max = max(TIK_L1_U_ADMM_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L1_U_mus,TIK_L1_U_ADMM_ISNRS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ISNR(\mu)  [dB]'); 
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TIK_L1_U_ADMM_ISSIMS); plot_y_max = max(TIK_L1_U_ADMM_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L1_U_mus,TIK_L1_U_ADMM_ISSIMS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ISSIM(\mu)');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TIK-L1-U-ADMM:    EFFICIENCY MEASURES  versus  \mu');
            yyaxis left;
            plot_y_min = min(TIK_L1_U_ADMM_ITRS); plot_y_max = max(TIK_L1_U_ADMM_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L1_U_mus,TIK_L1_U_ADMM_ITRS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ITERS(\mu)');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TIK_L1_U_ADMM_CPU_TIMES); plot_y_max = max(TIK_L1_U_ADMM_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L1_U_mus,TIK_L1_U_ADMM_CPU_TIMES,'-o','markersize',3);
            xlabel('\mu'); ylabel('CPU-TIME(\mu)  [secs]');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_TAU_STAR_VS_MU == 1 )
            figure('Name','Figure 34:  TIK-L1-U-ADMM:  TAU*  versus  MU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));            
            plot_y_min = min(TIK_L1_U_ADMM_TAU_STARS); plot_y_max = max(TIK_L1_U_ADMM_TAU_STARS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L1_U_mus,TIK_L1_U_ADMM_TAU_STARS,'-om','markersize',3);
            xlabel('\mu'); ylabel('\tau^*(\mu)');
            title('TIK-L1-U-ADMM:    TAU*  versus  \mu');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_PERFORMANCE_GRAPHS_VS_TAU_STAR == 1 )
            tau_stars_min = min(TIK_L1_U_ADMM_TAU_STARS); tau_stars_max = max(TIK_L1_U_ADMM_TAU_STARS);
            figure('Name','Figure 35:  TIK-L1-U-ADMM:  PERFORMANCE GRAPHS  versus  TAU*(MU)','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TIK-L1-U-ADMM:    ACCURACY MEASURES  versus  \tau^*(\mu)');
            yyaxis left;
            plot_y_min = min(TIK_L1_U_ADMM_ISNRS); plot_y_max = max(TIK_L1_U_ADMM_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L1_U_ADMM_TAU_STARS,TIK_L1_U_ADMM_ISNRS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ISNR(\tau^*(\mu))  [dB]'); 
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TIK_L1_U_ADMM_ISSIMS); plot_y_max = max(TIK_L1_U_ADMM_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L1_U_ADMM_TAU_STARS,TIK_L1_U_ADMM_ISSIMS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ISSIM(\tau^*(\mu))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TIK-L1-U-ADMM:    EFFICIENCY MEASURES  versus  \tau^*(\mu)');
            yyaxis left;
            plot_y_min = min(TIK_L1_U_ADMM_ITRS); plot_y_max = max(TIK_L1_U_ADMM_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L1_U_ADMM_TAU_STARS,TIK_L1_U_ADMM_ITRS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ITERS(\tau^*(\mu))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TIK_L1_U_ADMM_CPU_TIMES); plot_y_max = max(TIK_L1_U_ADMM_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L1_U_ADMM_TAU_STARS,TIK_L1_U_ADMM_CPU_TIMES,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('CPU-TIME(\tau^*(\mu))  [secs]');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
        end
    end
end

%--------------------------------------------------------------------------

if ( TEST_TV_L2_U_ADMM == 1 )
    % compute (in order to print it) tau_star of u0
    Au0_b       = Ax(u0) - b;
    u0_tau_star = norm( Au0_b(:) ) / (sqrt(d)*n_sigma);
    % print preliminaries of model/algorithm/initial iterate in the command window
    fprintf('\n\n\nIMAGE %s by TV-L2-U MODEL solved by ADMM:\n',INV_PROB_NAME_LONG);
    if ( TV_L2_U_mus_n == 1 ), itrs_debugs = 1; else, itrs_debugs = 0; 
        fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);    
    end
    % cycle over all the selected different mu values
    TV_L2_U_ADMM_BEST_ISNR  = -10^5; 
    TV_L2_U_ADMM_BEST_ISSIM = -10^5;
    for mu_i = 1:TV_L2_U_mus_n
        mu = TV_L2_U_mus(mu_i);
        % for the current value of the regularization parameter mu, 
        % compute the reconstructed image u_star(mu) by solving the model by ADMM        
        t0 = cputime;
        
        if ( TEST_REST_OR_INPT == 0 )
            [u_star,itrs] = REST_TV_L2_U_ADMM(b,blur_k,mu,TV_L2_U_ADMM_beta_t,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
        else
            if ( FAST_INPT == 1 )
                [u_star,itrs] = INPT_TV_L2_U_ADMM_FAST(b,M,mu,TV_L2_U_ADMM_beta_t,TV_L2_U_ADMM_beta_t,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
            else
                [u_star,itrs] = INPT_TV_L2_U_ADMM(b,M,mu,TV_L2_U_ADMM_beta_t,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
            end
        end
        
        % compute/store the current reconstruction scalar performance measures
        TV_L2_U_ADMM_CPU_TIMES(mu_i) = cputime - t0;
        TV_L2_U_ADMM_ISNRS(mu_i)     = compute_snr(u_star,u_true) - u0_SNR;
        TV_L2_U_ADMM_ISSIMS(mu_i)    = ssim(u_star,u_true) - u0_SSIM;
        TV_L2_U_ADMM_ITRS(mu_i)      = itrs;
        Au_star_b                    = Ax(u_star) - b;
        TV_L2_U_ADMM_TAU_STARS(mu_i) = norm( Au_star_b(:) ) / (sqrt(d)*n_sigma);
        % update the best reconstruction results in terms of accuracy measures ISNR and SSIM 
        if ( TV_L2_U_ADMM_ISNRS(mu_i) > TV_L2_U_ADMM_BEST_ISNR )
            TV_L2_U_ADMM_BEST_ISNR         = TV_L2_U_ADMM_ISNRS(mu_i);
            TV_L2_U_ADMM_BEST_ISNR_u_star  = u_star;
            TV_L2_U_ADMM_BEST_ISNR_mu_i    = mu_i;
        end
        if ( TV_L2_U_ADMM_ISSIMS(mu_i) > TV_L2_U_ADMM_BEST_ISSIM )
            TV_L2_U_ADMM_BEST_ISSIM        = TV_L2_U_ADMM_ISSIMS(mu_i);
            TV_L2_U_ADMM_BEST_ISSIM_u_star = u_star;
            TV_L2_U_ADMM_BEST_ISSIM_mu_i   = mu_i;
        end
        % print in the Matlab command window the current reconstruction results
        if ( TV_L2_U_mus_n == 1 )
            fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);
        end
        fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,  [ISNR,ISSIM] = [%6.3f , %6.4f]',...
            mu_i,TV_L2_U_mus_n,mu,itrs,TV_L2_U_ADMM_TAU_STARS(mu_i),TV_L2_U_ADMM_ISNRS(mu_i),TV_L2_U_ADMM_ISSIMS(mu_i));
        % show in a Matlab figure the current reconstruction results
        if ( SHOW_ALL_RESTORED_IMAGES == 1 )
            if ( mu_i == 1 )
                figure('Name','Figure 21:  TV-L2-U-ADMM:   mu-DEPENDENT RESULTS','NumberTitle','off');
                set(gcf,'Position',get(0,'ScreenSize'));
            end
            subplot(2,2,1)
            imshow(uint8(IMS_VIS_SF*u_true));
            title(sprintf('IM%02d (%d x %d):   u_{true}',im_id,h,w));
            subplot(2,2,3)
            if ( TEST_REST_OR_INPT == 0 )
                imshow(uint8(IMS_VIS_SF*b));
            else
                imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
            end
            title(sprintf('b = N(A u_{true}):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM));
            subplot(2,2,2)
            imshow(uint8(IMS_VIS_SF*u_star));
            title(sprintf('u^*(mu = %7.2f):  [ISNR,ISSIM] = [%6.3f,%6.4f]',mu,TV_L2_U_ADMM_ISNRS(mu_i),TV_L2_U_ADMM_ISSIMS(mu_i)));
            subplot(2,2,4)
            imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star) )); 
            title(sprintf('| u^*(mu = %7.2f) - u_{true} |',mu));
            pause(0.1);
        end
        % save on the hard disk the current reconstruction results
        if ( SAVE_ALL_RESTORED_IMAGES == 1 )
            imwrite(uint8(IMS_VIS_SF*u_star),sprintf('%s/%s_IM%02d_D2_TV_L2_U_ADMM_A_REST_mu_%07.2f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,mu));
            imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star)),sprintf('%s/%s_IM%02d_D2_TV_L2_U_ADMM_B_AERR_mu_%07.2f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,mu));
        end
    end
    % print in the Matlab command window the best reconstruction results in terms of ISNR and ISSIM
    TV_L2_U_ADMM_BEST_ISNR_mu        = TV_L2_U_mus(TV_L2_U_ADMM_BEST_ISNR_mu_i);
    TV_L2_U_ADMM_BEST_ISSIM_mu       = TV_L2_U_mus(TV_L2_U_ADMM_BEST_ISSIM_mu_i);
    TV_L2_U_ADMM_BEST_ISNR_itrs      = TV_L2_U_ADMM_ITRS(TV_L2_U_ADMM_BEST_ISNR_mu_i);
    TV_L2_U_ADMM_BEST_ISSIM_itrs     = TV_L2_U_ADMM_ITRS(TV_L2_U_ADMM_BEST_ISSIM_mu_i);
    TV_L2_U_ADMM_BEST_ISNR_tau_star  = TV_L2_U_ADMM_TAU_STARS(TV_L2_U_ADMM_BEST_ISNR_mu_i);
    TV_L2_U_ADMM_BEST_ISSIM_tau_star = TV_L2_U_ADMM_TAU_STARS(TV_L2_U_ADMM_BEST_ISSIM_mu_i);        
    fprintf('\nBEST ISNR and ISSIM RESULTS:');
    fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,   BEST ISNR   =  %6.3f',...
        TV_L2_U_ADMM_BEST_ISNR_mu_i,TV_L2_U_mus_n,TV_L2_U_ADMM_BEST_ISNR_mu,...
        TV_L2_U_ADMM_BEST_ISNR_itrs,TV_L2_U_ADMM_BEST_ISNR_tau_star,TV_L2_U_ADMM_BEST_ISNR);
    fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,   BEST ISSIM  =   %6.4f',...
        TV_L2_U_ADMM_BEST_ISSIM_mu_i,TV_L2_U_mus_n,TV_L2_U_ADMM_BEST_ISSIM_mu,...
        TV_L2_U_ADMM_BEST_ISSIM_itrs,TV_L2_U_ADMM_BEST_ISSIM_tau_star,TV_L2_U_ADMM_BEST_ISSIM);
    % show in a Matlab figure the best reconstruction results in terms of ISNR and ISSIM
    if ( SHOW_BEST_RESTORED_IMAGES == 1 )
        figure('Name','Figure 22:  TV-L2-U-ADMM:  BEST RESULTS in terms of ISNR and ISSIM','NumberTitle','off');
        set(gcf,'Position',get(0,'ScreenSize'));
        subplot(2,3,1)
        imshow(uint8(IMS_VIS_SF*u_true));
        title(sprintf('IM%02d (%d x %d):   u_{true}',im_id,h,w));
        subplot(2,3,4)
        if ( TEST_REST_OR_INPT == 0 )
            imshow(uint8(IMS_VIS_SF*b));
        else
            imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
        end
        title(sprintf('b = N(A u_{true}):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM));
        subplot(2,3,2)
        imshow(uint8(IMS_VIS_SF*TV_L2_U_ADMM_BEST_ISNR_u_star));
        title(sprintf('u^*(mu = %7.2f):  BEST ISNR = %6.3f',TV_L2_U_ADMM_BEST_ISNR_mu,TV_L2_U_ADMM_BEST_ISNR));
        subplot(2,3,5)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L2_U_ADMM_BEST_ISNR_u_star) )); 
        title(sprintf('| u^*(mu = %7.2f) - u_{true} |',TV_L2_U_ADMM_BEST_ISNR_mu));
        subplot(2,3,3)
        imshow(uint8(IMS_VIS_SF*TV_L2_U_ADMM_BEST_ISSIM_u_star));
        title(sprintf('u^*(mu = %7.2f):  BEST ISSIM = %6.4f',TV_L2_U_ADMM_BEST_ISSIM_mu,TV_L2_U_ADMM_BEST_ISSIM));
        subplot(2,3,6)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L2_U_ADMM_BEST_ISSIM_u_star) )); 
        title(sprintf('| u^*(mu = %7.2f) - u_{true} |',TV_L2_U_ADMM_BEST_ISSIM_mu));
    end
    % save on the hard disk the best reconstruction results in terms of ISNR and ISSIM
    if ( SAVE_BEST_RESTORED_IMAGES == 1 )
        imwrite(uint8(IMS_VIS_SF*TV_L2_U_ADMM_BEST_ISNR_u_star),sprintf('%s/%s_IM%02d_D2_TV_L2_U_ADMM_A_REST_mu_%07.2f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L2_U_ADMM_BEST_ISNR_mu));
        imwrite(uint8(IMS_VIS_SF*TV_L2_U_ADMM_BEST_ISSIM_u_star),sprintf('%s/%s_IM%02d_D2_TV_L2_U_ADMM_A_REST_mu_%07.2f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L2_U_ADMM_BEST_ISSIM_mu));
        imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L2_U_ADMM_BEST_ISNR_u_star)),sprintf('%s/%s_IM%02d_D2_TV_L2_U_ADMM_B_AERR_mu_%07.2f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L2_U_ADMM_BEST_ISNR_mu));
        imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L2_U_ADMM_BEST_ISSIM_u_star)),sprintf('%s/%s_IM%02d_D2_TV_L2_U_ADMM_B_AERR_mu_%07.2f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L2_U_ADMM_BEST_ISSIM_mu));
    end
    % show in a Matlab figure the performance graphs
    % show...only if we have selected at least two different mu values
    mus_min = min(TV_L2_U_mus); mus_max = max(TV_L2_U_mus);
    if ( mus_max > mus_min )
        if ( SHOW_PERFORMANCE_GRAPHS_VS_MU == 1 )
            figure('Name','Figure 23:  TV-L2-U-ADMM:  PERFORMANCE GRAPHS  versus  MU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TV-L2-U-ADMM:    ACCURACY MEASURES  versus  \mu');
            yyaxis left;
            plot_y_min = min(TV_L2_U_ADMM_ISNRS); plot_y_max = max(TV_L2_U_ADMM_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_U_mus,TV_L2_U_ADMM_ISNRS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ISNR(\mu)  [dB]'); 
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L2_U_ADMM_ISSIMS); plot_y_max = max(TV_L2_U_ADMM_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_U_mus,TV_L2_U_ADMM_ISSIMS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ISSIM(\mu)');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TV-L2-U-ADMM:    EFFICIENCY MEASURES  versus  \mu');
            yyaxis left;
            plot_y_min = min(TV_L2_U_ADMM_ITRS); plot_y_max = max(TV_L2_U_ADMM_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_U_mus,TV_L2_U_ADMM_ITRS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ITERS(\mu)');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L2_U_ADMM_CPU_TIMES); plot_y_max = max(TV_L2_U_ADMM_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_U_mus,TV_L2_U_ADMM_CPU_TIMES,'-o','markersize',3);
            xlabel('\mu'); ylabel('CPU-TIME(\mu)  [secs]');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_TAU_STAR_VS_MU == 1 )
            figure('Name','Figure 24:  TV-L2-U-ADMM:  TAU*  versus  MU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));            
            plot_y_min = min(TV_L2_U_ADMM_TAU_STARS); plot_y_max = max(TV_L2_U_ADMM_TAU_STARS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_U_mus,TV_L2_U_ADMM_TAU_STARS,'-om','markersize',3);
            xlabel('\mu'); ylabel('\tau^*(\mu)');
            title('TV-L2-U-ADMM:    TAU*  versus  \mu');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_PERFORMANCE_GRAPHS_VS_TAU_STAR == 1 )
            tau_stars_min = min(TV_L2_U_ADMM_TAU_STARS); tau_stars_max = max(TV_L2_U_ADMM_TAU_STARS);
            figure('Name','Figure 25:  TV-L2-U-ADMM:  PERFORMANCE GRAPHS  versus  TAU*(MU)','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TV-L2-U-ADMM:    ACCURACY MEASURES  versus  \tau^*(\mu)');
            yyaxis left;
            plot_y_min = min(TV_L2_U_ADMM_ISNRS); plot_y_max = max(TV_L2_U_ADMM_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_U_ADMM_TAU_STARS,TV_L2_U_ADMM_ISNRS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ISNR(\tau^*(\mu))  [dB]'); 
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L2_U_ADMM_ISSIMS); plot_y_max = max(TV_L2_U_ADMM_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_U_ADMM_TAU_STARS,TV_L2_U_ADMM_ISSIMS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ISSIM(\tau^*(\mu))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TV-L2-U-ADMM:    EFFICIENCY MEASURES  versus  \tau^*(\mu)');
            yyaxis left;
            plot_y_min = min(TV_L2_U_ADMM_ITRS); plot_y_max = max(TV_L2_U_ADMM_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_U_ADMM_TAU_STARS,TV_L2_U_ADMM_ITRS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ITERS(\tau^*(\mu))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L2_U_ADMM_CPU_TIMES); plot_y_max = max(TV_L2_U_ADMM_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_U_ADMM_TAU_STARS,TV_L2_U_ADMM_CPU_TIMES,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('CPU-TIME(\tau^*(\mu))  [secs]');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
        end
    end
end

%--------------------------------------------------------------------------

if ( TEST_TV_L1_U_ADMM == 1 )
    % compute (in order to print it) tau_star of u0
    Au0_b       = Ax(u0) - b;
    u0_tau_star = norm( Au0_b(:) ) / (sqrt(d)*n_sigma);
    % print preliminaries of model/algorithm/initial iterate in the command window
    fprintf('\n\n\nIMAGE %s by TV-L1-U MODEL solved by ADMM:\n',INV_PROB_NAME_LONG);
    if ( TV_L1_U_mus_n == 1 ), itrs_debugs = 1; else, itrs_debugs = 0; 
        fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);    
    end
    % cycle over all the selected different mu values
    TV_L1_U_ADMM_BEST_ISNR  = -10^5; 
    TV_L1_U_ADMM_BEST_ISSIM = -10^5;
    for mu_i = 1:TV_L1_U_mus_n
        mu = TV_L1_U_mus(mu_i);
        % for the current value of the regularization parameter mu, 
        % compute the reconstructed image u_star(mu) by solving the model by ADMM        
        t0 = cputime;
        
        if ( TEST_REST_OR_INPT == 0 )
            [u_star,itrs] = REST_TV_L1_U_ADMM(b,blur_k,mu,TV_L1_U_ADMM_beta_t,TV_L1_U_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
        else
            if ( FAST_INPT == 1 )
                [u_star,itrs] = INPT_TV_L1_U_ADMM_FAST(b,M,mu,TV_L1_U_ADMM_beta_t,TV_L1_U_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
            else
                [u_star,itrs] = INPT_TV_L1_U_ADMM(b,M,mu,TV_L1_U_ADMM_beta_t,TV_L1_U_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
            end
        end
        
        % compute/store the current reconstruction scalar performance measures
        TV_L1_U_ADMM_CPU_TIMES(mu_i) = cputime - t0;
        TV_L1_U_ADMM_ISNRS(mu_i)     = compute_snr(u_star,u_true) - u0_SNR;
        TV_L1_U_ADMM_ISSIMS(mu_i)    = ssim(u_star,u_true) - u0_SSIM;
        TV_L1_U_ADMM_ITRS(mu_i)      = itrs;
        Au_star_b                    = Ax(u_star) - b;
        TV_L1_U_ADMM_TAU_STARS(mu_i) = norm( Au_star_b(:) ) / (sqrt(d)*n_sigma);
        % update the best reconstruction results in terms of accuracy measures ISNR and SSIM 
        if ( TV_L1_U_ADMM_ISNRS(mu_i) > TV_L1_U_ADMM_BEST_ISNR )
            TV_L1_U_ADMM_BEST_ISNR         = TV_L1_U_ADMM_ISNRS(mu_i);
            TV_L1_U_ADMM_BEST_ISNR_u_star  = u_star;
            TV_L1_U_ADMM_BEST_ISNR_mu_i    = mu_i;
        end
        if ( TV_L1_U_ADMM_ISSIMS(mu_i) > TV_L1_U_ADMM_BEST_ISSIM )
            TV_L1_U_ADMM_BEST_ISSIM        = TV_L1_U_ADMM_ISSIMS(mu_i);
            TV_L1_U_ADMM_BEST_ISSIM_u_star = u_star;
            TV_L1_U_ADMM_BEST_ISSIM_mu_i   = mu_i;
        end
        % print in the Matlab command window the current reconstruction results
        if ( TV_L1_U_mus_n == 1 )
            fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);
        end
        fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,  [ISNR,ISSIM] = [%6.3f , %6.4f]',...
            mu_i,TV_L1_U_mus_n,mu,itrs,TV_L1_U_ADMM_TAU_STARS(mu_i),TV_L1_U_ADMM_ISNRS(mu_i),TV_L1_U_ADMM_ISSIMS(mu_i));
        % show in a Matlab figure the current reconstruction results
        if ( SHOW_ALL_RESTORED_IMAGES == 1 )
            if ( mu_i == 1 )
                figure('Name','Figure 41:  TV-L1-U-ADMM:   mu-DEPENDENT RESULTS','NumberTitle','off');
                set(gcf,'Position',get(0,'ScreenSize'));
            end
            subplot(2,2,1)
            imshow(uint8(IMS_VIS_SF*u_true));
            title(sprintf('IM%02d (%d x %d):   u_{true}',im_id,h,w));
            subplot(2,2,3)
            if ( TEST_REST_OR_INPT == 0 )
                imshow(uint8(IMS_VIS_SF*b));
            else
                imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
            end
            title(sprintf('b = N(A u_{true}):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM));
            subplot(2,2,2)
            imshow(uint8(IMS_VIS_SF*u_star));
            title(sprintf('u^*(mu = %7.2f):  [ISNR,ISSIM] = [%6.3f,%6.4f]',mu,TV_L1_U_ADMM_ISNRS(mu_i),TV_L1_U_ADMM_ISSIMS(mu_i)));
            subplot(2,2,4)
            imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star) )); 
            title(sprintf('| u^*(mu = %7.2f) - u_{true} |',mu));
            pause(0.1);
        end
        % save on the hard disk the current reconstruction results
        if ( SAVE_ALL_RESTORED_IMAGES == 1 )
            imwrite(uint8(IMS_VIS_SF*u_star),sprintf('%s/%s_IM%02d_D4_TV_L1_U_ADMM_A_REST_mu_%07.2f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,mu));
            imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star)),sprintf('%s/%s_IM%02d_D4_TV_L1_U_ADMM_B_AERR_mu_%07.2f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,mu));
        end
    end
    % print in the Matlab command window the best reconstruction results in terms of ISNR and ISSIM
    TV_L1_U_ADMM_BEST_ISNR_mu        = TV_L1_U_mus(TV_L1_U_ADMM_BEST_ISNR_mu_i);
    TV_L1_U_ADMM_BEST_ISSIM_mu       = TV_L1_U_mus(TV_L1_U_ADMM_BEST_ISSIM_mu_i);
    TV_L1_U_ADMM_BEST_ISNR_itrs      = TV_L1_U_ADMM_ITRS(TV_L1_U_ADMM_BEST_ISNR_mu_i);
    TV_L1_U_ADMM_BEST_ISSIM_itrs     = TV_L1_U_ADMM_ITRS(TV_L1_U_ADMM_BEST_ISSIM_mu_i);
    TV_L1_U_ADMM_BEST_ISNR_tau_star  = TV_L1_U_ADMM_TAU_STARS(TV_L1_U_ADMM_BEST_ISNR_mu_i);
    TV_L1_U_ADMM_BEST_ISSIM_tau_star = TV_L1_U_ADMM_TAU_STARS(TV_L1_U_ADMM_BEST_ISSIM_mu_i);        
    fprintf('\nBEST ISNR and ISSIM RESULTS:');
    fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,   BEST ISNR   =  %6.3f',...
        TV_L1_U_ADMM_BEST_ISNR_mu_i,TV_L1_U_mus_n,TV_L1_U_ADMM_BEST_ISNR_mu,...
        TV_L1_U_ADMM_BEST_ISNR_itrs,TV_L1_U_ADMM_BEST_ISNR_tau_star,TV_L1_U_ADMM_BEST_ISNR);
    fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,   BEST ISSIM  =   %6.4f',...
        TV_L1_U_ADMM_BEST_ISSIM_mu_i,TV_L1_U_mus_n,TV_L1_U_ADMM_BEST_ISSIM_mu,...
        TV_L1_U_ADMM_BEST_ISSIM_itrs,TV_L1_U_ADMM_BEST_ISSIM_tau_star,TV_L1_U_ADMM_BEST_ISSIM);
    % show in a Matlab figure the best reconstruction results in terms of ISNR and ISSIM
    if ( SHOW_BEST_RESTORED_IMAGES == 1 )
        figure('Name','Figure 42:  TV-L1-U-ADMM:  BEST RESULTS in terms of ISNR and ISSIM','NumberTitle','off');
        set(gcf,'Position',get(0,'ScreenSize'));
        subplot(2,3,1)
        imshow(uint8(IMS_VIS_SF*u_true));
        title(sprintf('IM%02d (%d x %d):   u_{true}',im_id,h,w));
        subplot(2,3,4)
        if ( TEST_REST_OR_INPT == 0 )
            imshow(uint8(IMS_VIS_SF*b));
        else
            imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
        end
        title(sprintf('b = N(A u_{true}):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM));
        subplot(2,3,2)
        imshow(uint8(IMS_VIS_SF*TV_L1_U_ADMM_BEST_ISNR_u_star));
        title(sprintf('u^*(mu = %7.2f):  BEST ISNR = %6.3f',TV_L1_U_ADMM_BEST_ISNR_mu,TV_L1_U_ADMM_BEST_ISNR));
        subplot(2,3,5)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L1_U_ADMM_BEST_ISNR_u_star) )); 
        title(sprintf('| u^*(mu = %7.2f) - u_{true} |',TV_L1_U_ADMM_BEST_ISNR_mu));
        subplot(2,3,3)
        imshow(uint8(IMS_VIS_SF*TV_L1_U_ADMM_BEST_ISSIM_u_star));
        title(sprintf('u^*(mu = %7.2f):  BEST ISSIM = %6.4f',TV_L1_U_ADMM_BEST_ISSIM_mu,TV_L1_U_ADMM_BEST_ISSIM));
        subplot(2,3,6)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L1_U_ADMM_BEST_ISSIM_u_star) )); 
        title(sprintf('| u^*(mu = %7.2f) - u_{true} |',TV_L1_U_ADMM_BEST_ISSIM_mu));
    end
    % save on the hard disk the best reconstruction results in terms of ISNR and ISSIM
    if ( SAVE_BEST_RESTORED_IMAGES == 1 )
        imwrite(uint8(IMS_VIS_SF*TV_L1_U_ADMM_BEST_ISNR_u_star),sprintf('%s/%s_IM%02d_D4_TV_L1_U_ADMM_A_REST_mu_%07.2f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L1_U_ADMM_BEST_ISNR_mu));
        imwrite(uint8(IMS_VIS_SF*TV_L1_U_ADMM_BEST_ISSIM_u_star),sprintf('%s/%s_IM%02d_D4_TV_L1_U_ADMM_A_REST_mu_%07.2f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L1_U_ADMM_BEST_ISSIM_mu));
        imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L1_U_ADMM_BEST_ISNR_u_star)),sprintf('%s/%s_IM%02d_D4_TV_L1_U_ADMM_B_AERR_mu_%07.2f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L1_U_ADMM_BEST_ISNR_mu));
        imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L1_U_ADMM_BEST_ISSIM_u_star)),sprintf('%s/%s_IM%02d_D4_TV_L1_U_ADMM_B_AERR_mu_%07.2f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L1_U_ADMM_BEST_ISSIM_mu));
    end
    % show in a Matlab figure the performance graphs
    % show...only if we have selected at least two different mu values
    mus_min = min(TV_L1_U_mus); mus_max = max(TV_L1_U_mus);
    if ( mus_max > mus_min )
        if ( SHOW_PERFORMANCE_GRAPHS_VS_MU == 1 )
            figure('Name','Figure 43:  TV-L1-U-ADMM:  PERFORMANCE GRAPHS  versus  MU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TV-L1-U-ADMM:    ACCURACY MEASURES  versus  \mu');
            yyaxis left;
            plot_y_min = min(TV_L1_U_ADMM_ISNRS); plot_y_max = max(TV_L1_U_ADMM_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L1_U_mus,TV_L1_U_ADMM_ISNRS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ISNR(\mu)  [dB]'); 
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L1_U_ADMM_ISSIMS); plot_y_max = max(TV_L1_U_ADMM_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L1_U_mus,TV_L1_U_ADMM_ISSIMS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ISSIM(\mu)');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TV-L1-U-ADMM:    EFFICIENCY MEASURES  versus  \mu');
            yyaxis left;
            plot_y_min = min(TV_L1_U_ADMM_ITRS); plot_y_max = max(TV_L1_U_ADMM_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L1_U_mus,TV_L1_U_ADMM_ITRS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ITERS(\mu)');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L1_U_ADMM_CPU_TIMES); plot_y_max = max(TV_L1_U_ADMM_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L1_U_mus,TV_L1_U_ADMM_CPU_TIMES,'-o','markersize',3);
            xlabel('\mu'); ylabel('CPU-TIME(\mu)  [secs]');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_TAU_STAR_VS_MU == 1 )
            figure('Name','Figure 44:  TV-L1-U-ADMM:  TAU*  versus  MU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));            
            plot_y_min = min(TV_L1_U_ADMM_TAU_STARS); plot_y_max = max(TV_L1_U_ADMM_TAU_STARS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L1_U_mus,TV_L1_U_ADMM_TAU_STARS,'-om','markersize',3);
            xlabel('\mu'); ylabel('\tau^*(\mu)');
            title('TV-L1-U-ADMM:    TAU*  versus  \mu');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_PERFORMANCE_GRAPHS_VS_TAU_STAR == 1 )
            tau_stars_min = min(TV_L1_U_ADMM_TAU_STARS); tau_stars_max = max(TV_L1_U_ADMM_TAU_STARS);
            figure('Name','Figure 45:  TV-L1-U-ADMM:  PERFORMANCE GRAPHS  versus  TAU*(MU)','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TV-L1-U-ADMM:    ACCURACY MEASURES  versus  \tau^*(\mu)');
            yyaxis left;
            plot_y_min = min(TV_L1_U_ADMM_ISNRS); plot_y_max = max(TV_L1_U_ADMM_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L1_U_ADMM_TAU_STARS,TV_L1_U_ADMM_ISNRS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ISNR(\tau^*(\mu))  [dB]'); 
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L1_U_ADMM_ISSIMS); plot_y_max = max(TV_L1_U_ADMM_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L1_U_ADMM_TAU_STARS,TV_L1_U_ADMM_ISSIMS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ISSIM(\tau^*(\mu))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TV-L1-U-ADMM:    EFFICIENCY MEASURES  versus  \tau^*(\mu)');
            yyaxis left;
            plot_y_min = min(TV_L1_U_ADMM_ITRS); plot_y_max = max(TV_L1_U_ADMM_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L1_U_ADMM_TAU_STARS,TV_L1_U_ADMM_ITRS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ITERS(\tau^*(\mu))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L1_U_ADMM_CPU_TIMES); plot_y_max = max(TV_L1_U_ADMM_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L1_U_ADMM_TAU_STARS,TV_L1_U_ADMM_CPU_TIMES,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('CPU-TIME(\tau^*(\mu))  [secs]');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
        end
    end
end

%--------------------------------------------------------------------------

if ( TEST_TIK_L15_U_ADMM == 1 )
    profile on
    % compute (in order to print it) tau_star of u0
    Au0_b       = Ax(u0) - b;
    u0_tau_star = norm( Au0_b(:) ) / (sqrt(d)*n_sigma);
    % print preliminaries of model/algorithm/initial iterate in the command window
    fprintf('\n\n\nIMAGE %s by TIK-L15-U MODEL solved by ADMM:\n',INV_PROB_NAME_LONG);
    if ( TIK_L15_U_mus_n == 1 ), itrs_debugs = 1; else, itrs_debugs = 0; 
        fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);    
    end
    % cycle over all the selected different mu values
    TIK_L15_U_ADMM_BEST_ISNR  = -10^5; 
    TIK_L15_U_ADMM_BEST_ISSIM = -10^5;
    ITS = 0;
    for mu_i = 1:TIK_L15_U_mus_n
        mu = TIK_L15_U_mus(mu_i);
        % for the current value of the regularization parameter mu,  
        % compute the reconstructed image u_star(mu) by solving the model by ADMM        
        t0 = cputime;
        
        if ( TEST_REST_OR_INPT == 0 )
            [u_star,itrs] = REST_TIK_L15_U_ADMM(b,blur_k,mu,TIK_L15_U_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
        else
            if ( FAST_INPT == 1 )
                [u_star,itrs] = INPT_TIK_L15_U_ADMM_FAST(b,M,mu,TIK_L15_U_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th);
            else
                [u_star,itrs] = INPT_TIK_L15_U_ADMM(b,M,mu,TIK_L15_U_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
            end
        end

        ITS = ITS + itrs;

        % compute/store the current reconstruction scalar performance measures
        TIK_L15_U_ADMM_CPU_TIMES(mu_i) = cputime - t0;
        TIK_L15_U_ADMM_ISNRS(mu_i)     = compute_snr(u_star,u_true) - u0_SNR;
        TIK_L15_U_ADMM_ISSIMS(mu_i)    = ssim(u_star,u_true) - u0_SSIM;
        TIK_L15_U_ADMM_ITRS(mu_i)      = itrs;
        Au_star_b                     = Ax(u_star) - b;
        TIK_L15_U_ADMM_TAU_STARS(mu_i) = norm( Au_star_b(:) ) / (sqrt(d)*n_sigma);
        % update the best reconstruction results in terms of accuracy measures ISNR and SSIM 
        if ( TIK_L15_U_ADMM_ISNRS(mu_i) > TIK_L15_U_ADMM_BEST_ISNR )
            TIK_L15_U_ADMM_BEST_ISNR         = TIK_L15_U_ADMM_ISNRS(mu_i);
            TIK_L15_U_ADMM_BEST_ISNR_u_star  = u_star;
            TIK_L15_U_ADMM_BEST_ISNR_mu_i    = mu_i;
        end
        if ( TIK_L15_U_ADMM_ISSIMS(mu_i) > TIK_L15_U_ADMM_BEST_ISSIM )
            TIK_L15_U_ADMM_BEST_ISSIM        = TIK_L15_U_ADMM_ISSIMS(mu_i);
            TIK_L15_U_ADMM_BEST_ISSIM_u_star = u_star;
            TIK_L15_U_ADMM_BEST_ISSIM_mu_i   = mu_i;
        end
        % print in the Matlab command window the current reconstruction results
        if ( TIK_L15_U_mus_n == 1 )
            fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);
        end
        fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,  [ISNR,ISSIM] = [%6.3f , %6.4f]',...
            mu_i,TIK_L15_U_mus_n,mu,itrs,TIK_L15_U_ADMM_TAU_STARS(mu_i),TIK_L15_U_ADMM_ISNRS(mu_i),TIK_L15_U_ADMM_ISSIMS(mu_i));
        % show in a Matlab figure the current reconstruction results
        if ( SHOW_ALL_RESTORED_IMAGES == 1 )
            if ( mu_i == 1 )
                figure('Name','Figure 31:  TIK-L1-U-ADMM:   mu-DEPENDENT RESULTS','NumberTitle','off');
                set(gcf,'Position',get(0,'ScreenSize'));
            end
            subplot(2,2,1)
            imshow(uint8(IMS_VIS_SF*u_true));
            title(sprintf('IM%02d (%d x %d):   u_{true}',im_id,h,w));
            subplot(2,2,3)
            if ( TEST_REST_OR_INPT == 0 )
                imshow(uint8(IMS_VIS_SF*b));
            else
                imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
            end
            title(sprintf('b = N(A u_{true}):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM));
            subplot(2,2,2)
            imshow(uint8(IMS_VIS_SF*u_star));
            title(sprintf('u^*(mu = %7.2f):  [ISNR,ISSIM] = [%6.3f,%6.4f]',mu,TIK_L15_U_ADMM_ISNRS(mu_i),TIK_L15_U_ADMM_ISSIMS(mu_i)));
            subplot(2,2,4)
            imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star) )); 
            title(sprintf('| u^*(mu = %7.2f) - u_{true} |',mu));
            pause(0.1);
        end
        % save on the hard disk the current reconstruction results
        if ( SAVE_ALL_RESTORED_IMAGES == 1 )
            imwrite(uint8(IMS_VIS_SF*u_star),sprintf('%s/%s_IM%02d_D3_TIK_L15_U_ADMM_A_REST_mu_%07.2f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,mu));
            imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star)),sprintf('%s/%s_IM%02d_D3_TIK_L15_U_ADMM_B_AERR_mu_%07.2f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,mu));
        end
    end
    % print in the Matlab command window the best reconstruction results in terms of ISNR and ISSIM
    TIK_L15_U_ADMM_BEST_ISNR_mu        = TIK_L15_U_mus(TIK_L15_U_ADMM_BEST_ISNR_mu_i);
    TIK_L15_U_ADMM_BEST_ISSIM_mu       = TIK_L15_U_mus(TIK_L15_U_ADMM_BEST_ISSIM_mu_i);
    TIK_L15_U_ADMM_BEST_ISNR_itrs      = TIK_L15_U_ADMM_ITRS(TIK_L15_U_ADMM_BEST_ISNR_mu_i);
    TIK_L15_U_ADMM_BEST_ISSIM_itrs     = TIK_L15_U_ADMM_ITRS(TIK_L15_U_ADMM_BEST_ISSIM_mu_i);
    TIK_L15_U_ADMM_BEST_ISNR_tau_star  = TIK_L15_U_ADMM_TAU_STARS(TIK_L15_U_ADMM_BEST_ISNR_mu_i);
    TIK_L15_U_ADMM_BEST_ISSIM_tau_star = TIK_L15_U_ADMM_TAU_STARS(TIK_L15_U_ADMM_BEST_ISSIM_mu_i);  

    % Best ISNR and ISSIM results
      fprintf('\nBEST ISNR and ISSIM RESULTS:');
      fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,   BEST ISNR   =  %6.3f',...
          TIK_L15_U_ADMM_BEST_ISNR_mu_i,TIK_L15_U_mus_n,TIK_L15_U_ADMM_BEST_ISNR_mu,...
          TIK_L15_U_ADMM_BEST_ISNR_itrs,TIK_L15_U_ADMM_BEST_ISNR_tau_star,TIK_L15_U_ADMM_BEST_ISNR);
      fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,   BEST ISSIM  =   %6.4f',...
          TIK_L15_U_ADMM_BEST_ISSIM_mu_i,TIK_L15_U_mus_n,TIK_L15_U_ADMM_BEST_ISSIM_mu,...
          TIK_L15_U_ADMM_BEST_ISSIM_itrs,TIK_L15_U_ADMM_BEST_ISSIM_tau_star,TIK_L15_U_ADMM_BEST_ISSIM);

    % Show in a Matlab figure the best reconstruction results in terms of
    % ISNR and ISSIM:
    if ( SHOW_BEST_RESTORED_IMAGES == 1 )
        figure('Name','Figure 32:  TIK-L1-U-ADMM:  BEST RESULTS in terms of ISNR and ISSIM','NumberTitle','off');
        set(gcf,'Position',get(0,'ScreenSize'));
        subplot(2,3,1)
        imshow(uint8(IMS_VIS_SF*u_true));
        title(sprintf('IM%02d (%d x %d):   $u_{true}$',im_id,h,w),'Interpreter','Latex');
        subplot(2,3,4)
        if ( TEST_REST_OR_INPT == 0 )
            imshow(uint8(IMS_VIS_SF*b));
        else
            imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
        end
        title(sprintf('b = N(A $u_{true}$):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM),'Interpreter','Latex');
        subplot(2,3,2)
        imshow(uint8(IMS_VIS_SF*TIK_L15_U_ADMM_BEST_ISNR_u_star));
        title(sprintf('$u^{*}$ ($\\mu$ = %7.2f ):  BEST ISNR = %6.3f',TIK_L15_U_ADMM_BEST_ISNR_mu,TIK_L15_U_ADMM_BEST_ISNR),'Interpreter','Latex');
        subplot(2,3,5)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TIK_L15_U_ADMM_BEST_ISNR_u_star) )); 
        title(sprintf('$|$ $u^{*}$ ($\\mu$ = %7.2f ) - $u_{true}$ $|$',TIK_L15_U_ADMM_BEST_ISNR_mu),'interpreter','latex');
        subplot(2,3,3)
        imshow(uint8(IMS_VIS_SF*TIK_L15_U_ADMM_BEST_ISSIM_u_star));
        title(sprintf('$u^{*}(\\mu$ = %7.2f):  BEST ISSIM = %6.4f',TIK_L15_U_ADMM_BEST_ISSIM_mu,TIK_L15_U_ADMM_BEST_ISSIM),'Interpreter','Latex');
        subplot(2,3,6)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TIK_L15_U_ADMM_BEST_ISSIM_u_star) )); 
        title(sprintf('$|$ $u^{*}$($\\mu$ = %7.2f) - $u_{true}$ $|$',TIK_L15_U_ADMM_BEST_ISSIM_mu),'interpreter','Latex');
    end
    % save on the hard disk the best reconstruction results in terms of ISNR and ISSIM
    % if ( SAVE_BEST_RESTORED_IMAGES == 1 )
    %     imwrite(uint8(IMS_VIS_SF*TIK_L15_U_ADMM_BEST_ISNR_u_star),sprintf('%s/%s_IM%02d_D3_TIK_L15_U_ADMM_A_REST_mu_%07.2f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TIK_L15_U_ADMM_BEST_ISNR_mu));
    %     imwrite(uint8(IMS_VIS_SF*TIK_L15_U_ADMM_BEST_ISSIM_u_star),sprintf('%s/%s_IM%02d_D3_TIK_L15_U_ADMM_A_REST_mu_%07.2f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TIK_L15_U_ADMM_BEST_ISSIM_mu));
    %     imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TIK_L15_U_ADMM_BEST_ISNR_u_star)),sprintf('%s/%s_IM%02d_D3_TIK_L15_U_ADMM_B_AERR_mu_%07.2f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TIK_L15_U_ADMM_BEST_ISNR_mu));
    %     imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TIK_L15_U_ADMM_BEST_ISSIM_u_star)),sprintf('%s/%s_IM%02d_D3_TIK_L15_U_ADMM_B_AERR_mu_%07.2f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TIK_L15_U_ADMM_BEST_ISSIM_mu));
    % end

    % show in a Matlab figure the performance graphs 
    % only if we have selected at least two different mu values

    mus_min = min(TIK_L15_U_mus); mus_max = max(TIK_L15_U_mus);
    if ( mus_max > mus_min )
        if ( SHOW_PERFORMANCE_GRAPHS_VS_MU == 1 )
            figure('Name','Figure 33:  TIK-L_{1.5}-U-ADMM:  PERFORMANCE GRAPHS  versus  MU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('$TIK-L_{1.5}$-U-ADMM:    ACCURACY MEASURES  versus  $\mu$','Interpreter','Latex');
            yyaxis left;
            plot_y_min = min(TIK_L15_U_ADMM_ISNRS); plot_y_max = max(TIK_L15_U_ADMM_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L15_U_mus,TIK_L15_U_ADMM_ISNRS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ISNR(\mu)  [dB]'); 
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TIK_L15_U_ADMM_ISSIMS); plot_y_max = max(TIK_L15_U_ADMM_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L15_U_mus,TIK_L15_U_ADMM_ISSIMS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ISSIM(\mu)');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('$TIK-L_{1.5}$-U-ADMM:    EFFICIENCY MEASURES  versus  $\mu$','Interpreter','Latex');
            yyaxis left;
            plot_y_min = min(TIK_L15_U_ADMM_ITRS); plot_y_max = max(TIK_L15_U_ADMM_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L15_U_mus,TIK_L15_U_ADMM_ITRS,'-o','markersize',3);
            xlabel('\mu'); ylabel('ITERS(\mu)');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TIK_L15_U_ADMM_CPU_TIMES); plot_y_max = max(TIK_L15_U_ADMM_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L15_U_mus,TIK_L15_U_ADMM_CPU_TIMES,'-o','markersize',3);
            xlabel('\mu'); ylabel('CPU-TIME(\mu)  [secs]');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_TAU_STAR_VS_MU == 1 )
            figure('Name','Figure 34:  TIK-L_{1.5}-U-ADMM:  TAU*  versus  MU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));            
            plot_y_min = min(TIK_L15_U_ADMM_TAU_STARS); plot_y_max = max(TIK_L15_U_ADMM_TAU_STARS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L15_U_mus,TIK_L15_U_ADMM_TAU_STARS,'-om','markersize',3);
            xlabel('\mu'); ylabel('\tau^*(\mu)');
            title('TIK-L{1.5}-U-ADMM:    TAU*  versus  \mu','Interpreter','Latex');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_PERFORMANCE_GRAPHS_VS_TAU_STAR == 1 )
            tau_stars_min = min(TIK_L15_U_ADMM_TAU_STARS); tau_stars_max = max(TIK_L15_U_ADMM_TAU_STARS);
            figure('Name','Figure 35:  TIK-L_{1.5}-U-ADMM:  PERFORMANCE GRAPHS  versus  TAU*(MU)','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TIK-L_{1.5}-U-ADMM:    ACCURACY MEASURES  versus  \tau^*(\mu)','Interpreter','Latex');
            yyaxis left;
            plot_y_min = min(TIK_L15_U_ADMM_ISNRS); plot_y_max = max(TIK_L15_U_ADMM_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L15_U_ADMM_TAU_STARS,TIK_L15_U_ADMM_ISNRS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ISNR(\tau^*(\mu))  [dB]'); 
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TIK_L15_U_ADMM_ISSIMS); plot_y_max = max(TIK_L15_U_ADMM_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L15_U_ADMM_TAU_STARS,TIK_L15_U_ADMM_ISSIMS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ISSIM(\tau^*(\mu))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TIK-L_{1.5}-U-ADMM:    EFFICIENCY MEASURES  versus  \tau^*(\mu)');
            yyaxis left;
            plot_y_min = min(TIK_L15_U_ADMM_ITRS); plot_y_max = max(TIK_L15_U_ADMM_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L15_U_ADMM_TAU_STARS,TIK_L15_U_ADMM_ITRS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ITERS(\tau^*(\mu))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TIK_L15_U_ADMM_CPU_TIMES); plot_y_max = max(TIK_L15_U_ADMM_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TIK_L15_U_ADMM_TAU_STARS,TIK_L15_U_ADMM_CPU_TIMES,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('CPU-TIME(\tau^*(\mu))  [secs]');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
        end
    end
profile off
end

%-------------------------------------------------------------------------%

if ( TEST_TV_L15_U_ADMM == 1 )
    % compute (in order to print it) tau_star of u0
    Au0_b       = Ax(u0) - b;
    u0_tau_star = norm( Au0_b(:) ) / (sqrt(d)*n_sigma);
    % print preliminaries of model/algorithm/initial iterate in the command window
    fprintf('\n\n\nIMAGE %s by TV-L15-U MODEL solved by ADMM:\n',INV_PROB_NAME_LONG);
    if ( TV_L15_U_mus_n == 1 ), itrs_debugs = 1; else, itrs_debugs = 0; 
        fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);    
    end
    % cycle over all the selected different mu values
    TV_L15_U_ADMM_BEST_ISNR  = -10^5; 
    TV_L15_U_ADMM_BEST_ISSIM = -10^5;
    for mu_i = 1:TV_L15_U_mus_n
        mu = TV_L15_U_mus(mu_i);
        % for the current value of the regularization parameter mu, 
        % compute the reconstructed image u_star(mu) by solving the model by ADMM        
        t0 = cputime;
        
        if ( TEST_REST_OR_INPT == 0 )
            [u_star,itrs] = REST_TV_L15_U_ADMM(b,blur_k,mu,TV_L15_U_ADMM_beta_t,TV_L15_U_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
        else
            if ( FAST_INPT == 1 )
                [u_star,itrs] = INPT_TV_L15_U_ADMM_FAST(b,M,mu,TV_L15_U_ADMM_beta_t,TV_L15_U_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th);
            else
                [u_star,itrs] = INPT_TV_L15_U_ADMM(b,M,mu,TV_L15_U_ADMM_beta_t,TV_L15_U_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,n_sigma,u_true);
            end
        end
        
        % compute/store the current reconstruction scalar performance measures
        TV_L15_U_ADMM_CPU_TIMES(mu_i) = cputime - t0;
        TV_L15_U_ADMM_ISNRS(mu_i)     = compute_snr(u_star,u_true) - u0_SNR;
        TV_L15_U_ADMM_ISSIMS(mu_i)    = ssim(u_star,u_true) - u0_SSIM;
        TV_L15_U_ADMM_ITRS(mu_i)      = itrs;
        Au_star_b                    = Ax(u_star) - b;
        TV_L15_U_ADMM_TAU_STARS(mu_i) = norm( Au_star_b(:) ) / (sqrt(d)*n_sigma);
        % update the best reconstruction results in terms of accuracy measures ISNR and SSIM 
        if ( TV_L15_U_ADMM_ISNRS(mu_i) > TV_L15_U_ADMM_BEST_ISNR )
            TV_L15_U_ADMM_BEST_ISNR         = TV_L15_U_ADMM_ISNRS(mu_i);
            TV_L15_U_ADMM_BEST_ISNR_u_star  = u_star;
            TV_L15_U_ADMM_BEST_ISNR_mu_i    = mu_i;
        end
        if ( TV_L15_U_ADMM_ISSIMS(mu_i) > TV_L15_U_ADMM_BEST_ISSIM )
            TV_L15_U_ADMM_BEST_ISSIM        = TV_L15_U_ADMM_ISSIMS(mu_i);
            TV_L15_U_ADMM_BEST_ISSIM_u_star = u_star;
            TV_L15_U_ADMM_BEST_ISSIM_mu_i   = mu_i;
        end
        % print in the Matlab command window the current reconstruction results
        if ( TV_L15_U_mus_n == 1 )
            fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);
        end
        fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,  [ISNR,ISSIM] = [%6.3f , %6.4f]',...
            mu_i,TV_L15_U_mus_n,mu,itrs,TV_L15_U_ADMM_TAU_STARS(mu_i),TV_L15_U_ADMM_ISNRS(mu_i),TV_L15_U_ADMM_ISSIMS(mu_i));
        % show in a Matlab figure the current reconstruction results
        if ( SHOW_ALL_RESTORED_IMAGES == 1 )
            if ( mu_i == 1 )
                figure('Name','Figure 41:  TV-L1-U-ADMM:   mu-DEPENDENT RESULTS','NumberTitle','off');
                set(gcf,'Position',get(0,'ScreenSize'));
            end
            subplot(2,2,1)
            imshow(uint8(IMS_VIS_SF*u_true));
            title(sprintf('IM%02d (%d x %d):   u_{true}',im_id,h,w));
            subplot(2,2,3)
            if ( TEST_REST_OR_INPT == 0 )
                imshow(uint8(IMS_VIS_SF*b));
            else
                imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
            end
            title(sprintf('b = N(A u_{true}):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM));
            subplot(2,2,2)
            imshow(uint8(IMS_VIS_SF*u_star));
            title(sprintf('u^*(mu = %7.2f):  [ISNR,ISSIM] = [%6.3f,%6.4f]',mu,TV_L15_U_ADMM_ISNRS(mu_i),TV_L15_U_ADMM_ISSIMS(mu_i)));
            subplot(2,2,4)
            imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star) )); 
            title(sprintf('| u^*(mu = %7.2f) - u_{true} |',mu));
            pause(0.1);
        end
        % save on the hard disk the current reconstruction results
        if ( SAVE_ALL_RESTORED_IMAGES == 1 )
            imwrite(uint8(IMS_VIS_SF*u_star),sprintf('%s/%s_IM%02d_D4_TV_L15_U_ADMM_A_REST_mu_%07.2f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,mu));
            imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star)),sprintf('%s/%s_IM%02d_D4_TV_L15_U_ADMM_B_AERR_mu_%07.2f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,mu));
        end
    end
    % print in the Matlab command window the best reconstruction results in terms of ISNR and ISSIM
    TV_L15_U_ADMM_BEST_ISNR_mu        = TV_L15_U_mus(TV_L15_U_ADMM_BEST_ISNR_mu_i);
    TV_L15_U_ADMM_BEST_ISSIM_mu       = TV_L15_U_mus(TV_L15_U_ADMM_BEST_ISSIM_mu_i);
    TV_L15_U_ADMM_BEST_ISNR_itrs      = TV_L15_U_ADMM_ITRS(TV_L15_U_ADMM_BEST_ISNR_mu_i);
    TV_L15_U_ADMM_BEST_ISSIM_itrs     = TV_L15_U_ADMM_ITRS(TV_L15_U_ADMM_BEST_ISSIM_mu_i);
    TV_L15_U_ADMM_BEST_ISNR_tau_star  = TV_L15_U_ADMM_TAU_STARS(TV_L15_U_ADMM_BEST_ISNR_mu_i);
    TV_L15_U_ADMM_BEST_ISSIM_tau_star = TV_L15_U_ADMM_TAU_STARS(TV_L15_U_ADMM_BEST_ISSIM_mu_i);        
    fprintf('\nBEST ISNR and ISSIM RESULTS:');
    fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,   BEST ISNR   =  %6.3f',...
        TV_L15_U_ADMM_BEST_ISNR_mu_i,TV_L15_U_mus_n,TV_L15_U_ADMM_BEST_ISNR_mu,...
        TV_L15_U_ADMM_BEST_ISNR_itrs,TV_L15_U_ADMM_BEST_ISNR_tau_star,TV_L15_U_ADMM_BEST_ISNR);
    fprintf('\nmu(%02d/%02d) = %7.2f:  ITS = %04d,  tau* = %6.4f,   BEST ISSIM  =   %6.4f',...
        TV_L15_U_ADMM_BEST_ISSIM_mu_i,TV_L15_U_mus_n,TV_L15_U_ADMM_BEST_ISSIM_mu,...
        TV_L15_U_ADMM_BEST_ISSIM_itrs,TV_L15_U_ADMM_BEST_ISSIM_tau_star,TV_L15_U_ADMM_BEST_ISSIM);
    % show in a Matlab figure the best reconstruction results in terms of ISNR and ISSIM
    if ( SHOW_BEST_RESTORED_IMAGES == 1 )
        figure('Name','Figure 42:  TV-L_{1.5}-U-ADMM:  BEST RESULTS in terms of ISNR and ISSIM','NumberTitle','off');
        set(gcf,'Position',get(0,'ScreenSize'));
        subplot(2,3,1)
        imshow(uint8(IMS_VIS_SF*u_true));
        title(sprintf('IM%02d (%d x %d):   $u_{true}$',im_id,h,w),'interpreter','latex');
        subplot(2,3,4)
        if ( TEST_REST_OR_INPT == 0 )
            imshow(uint8(IMS_VIS_SF*b));
        else
            imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
        end
        title(sprintf('b = N(A $u_{true}$):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM),'Interpreter','latex');
        subplot(2,3,2)
        imshow(uint8(IMS_VIS_SF*TV_L15_U_ADMM_BEST_ISNR_u_star));
        title(sprintf('$u^{*}$($\\mu = %7.2f$):  BEST ISNR = %6.3f',TV_L15_U_ADMM_BEST_ISNR_mu,TV_L15_U_ADMM_BEST_ISNR),'Interpreter','latex');
        subplot(2,3,5)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L15_U_ADMM_BEST_ISNR_u_star) )); 
        title(sprintf('$|$ $u^{*}$($\\mu = %7.2f$) - $u_{true}$ $|$',TV_L15_U_ADMM_BEST_ISNR_mu),'Interpreter','latex');
        subplot(2,3,3)
        imshow(uint8(IMS_VIS_SF*TV_L15_U_ADMM_BEST_ISSIM_u_star));
        title(sprintf('$u^{*}$($\\mu = %7.2f$):  BEST ISSIM = %6.4f',TV_L15_U_ADMM_BEST_ISSIM_mu,TV_L15_U_ADMM_BEST_ISSIM),'Interpreter','latex');
        subplot(2,3,6)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L15_U_ADMM_BEST_ISSIM_u_star) )); 
        title(sprintf('$|$ $u^{*}$($\\mu = %7.2f$) - $u_{true}$ $|$',TV_L15_U_ADMM_BEST_ISSIM_mu),'Interpreter','latex');
    end
    % save on the hard disk the best reconstruction results in terms of ISNR and ISSIM
    if ( SAVE_BEST_RESTORED_IMAGES == 1 )
        imwrite(uint8(IMS_VIS_SF*TV_L15_U_ADMM_BEST_ISNR_u_star),sprintf('%s/%s_IM%02d_D4_TV_L15_U_ADMM_A_REST_mu_%07.2f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L15_U_ADMM_BEST_ISNR_mu));
        imwrite(uint8(IMS_VIS_SF*TV_L15_U_ADMM_BEST_ISSIM_u_star),sprintf('%s/%s_IM%02d_D4_TV_L15_U_ADMM_A_REST_mu_%07.2f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L15_U_ADMM_BEST_ISSIM_mu));
        imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L15_U_ADMM_BEST_ISNR_u_star)),sprintf('%s/%s_IM%02d_D4_TV_L15_U_ADMM_B_AERR_mu_%07.2f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L15_U_ADMM_BEST_ISNR_mu));
        imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L15_U_ADMM_BEST_ISSIM_u_star)),sprintf('%s/%s_IM%02d_D4_TV_L15_U_ADMM_B_AERR_mu_%07.2f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L15_U_ADMM_BEST_ISSIM_mu));
    end
    % show in a Matlab figure the performance graphs
    % show...only if we have selected at least two different mu values
    mus_min = min(TV_L15_U_mus); mus_max = max(TV_L15_U_mus);
    if ( mus_max > mus_min )
        if ( SHOW_PERFORMANCE_GRAPHS_VS_MU == 1 )
            figure('Name','Figure 43:  TV-L_{1.5}-U-ADMM:  PERFORMANCE GRAPHS  versus  MU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TV-$L_{1.5}$-U-ADMM:    ACCURACY MEASURES  versus  $\mu$','Interpreter','latex');
            yyaxis left;
            plot_y_min = min(TV_L15_U_ADMM_ISNRS); plot_y_max = max(TV_L15_U_ADMM_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L15_U_mus,TV_L15_U_ADMM_ISNRS,'-o','markersize',3);
            xlabel('$\mu$','Interpreter','latex'); ylabel('ISNR($\mu$)  [dB]','Interpreter','latex'); 
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L15_U_ADMM_ISSIMS); plot_y_max = max(TV_L15_U_ADMM_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L15_U_mus,TV_L15_U_ADMM_ISSIMS,'-o','markersize',3);
            xlabel('$\mu$','interpreter','latex'); ylabel('ISSIM($\mu$)','interpreter','latex');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TV-$L_{1.5}$-U-ADMM:    EFFICIENCY MEASURES  versus  $\mu$','interpreter','latex');
            yyaxis left;
            plot_y_min = min(TV_L15_U_ADMM_ITRS); plot_y_max = max(TV_L15_U_ADMM_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L15_U_mus,TV_L15_U_ADMM_ITRS,'-o','markersize',3);
            xlabel('$\mu$','interpreter','latex'); ylabel('ITERS($\mu$)','interpreter','latex');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L15_U_ADMM_CPU_TIMES); plot_y_max = max(TV_L15_U_ADMM_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L15_U_mus,TV_L15_U_ADMM_CPU_TIMES,'-o','markersize',3);
            xlabel('$\mu$','interpreter','latex'); ylabel('CPU-TIME($\mu$)  [secs]','Interpreter','latex');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_TAU_STAR_VS_MU == 1 )
            figure('Name','Figure 44:  TV-L1-U-ADMM:  TAU*  versus  MU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));            
            plot_y_min = min(TV_L15_U_ADMM_TAU_STARS); plot_y_max = max(TV_L15_U_ADMM_TAU_STARS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L15_U_mus,TV_L15_U_ADMM_TAU_STARS,'-om','markersize',3);
            xlabel('\mu'); ylabel('\tau^*(\mu)');
            title('TV-L1-U-ADMM:    TAU*  versus  \mu');
            axis([mus_min,mus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_PERFORMANCE_GRAPHS_VS_TAU_STAR == 1 )
            tau_stars_min = min(TV_L15_U_ADMM_TAU_STARS); tau_stars_max = max(TV_L15_U_ADMM_TAU_STARS);
            figure('Name','Figure 45:  TV-L1-U-ADMM:  PERFORMANCE GRAPHS  versus  TAU*(MU)','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TV-L1-U-ADMM:    ACCURACY MEASURES  versus  \tau^*(\mu)');
            yyaxis left;
            plot_y_min = min(TV_L15_U_ADMM_ISNRS); plot_y_max = max(TV_L15_U_ADMM_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L15_U_ADMM_TAU_STARS,TV_L15_U_ADMM_ISNRS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ISNR(\tau^*(\mu))  [dB]'); 
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L15_U_ADMM_ISSIMS); plot_y_max = max(TV_L15_U_ADMM_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L15_U_ADMM_TAU_STARS,TV_L15_U_ADMM_ISSIMS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ISSIM(\tau^*(\mu))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TV-L1-U-ADMM:    EFFICIENCY MEASURES  versus  \tau^*(\mu)');
            yyaxis left;
            plot_y_min = min(TV_L15_U_ADMM_ITRS); plot_y_max = max(TV_L15_U_ADMM_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L15_U_ADMM_TAU_STARS,TV_L15_U_ADMM_ITRS,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('ITERS(\tau^*(\mu))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L15_U_ADMM_CPU_TIMES); plot_y_max = max(TV_L15_U_ADMM_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L15_U_ADMM_TAU_STARS,TV_L15_U_ADMM_CPU_TIMES,'-o','markersize',3);
            xlabel('\tau^*(\mu)'); ylabel('CPU-TIME(\tau^*(\mu))  [secs]');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
        end
    end
end

%--------------------------------------------------------------------------

if ( TEST_TV_L2_DC_ADMM == 1 )
    % compute (in order to print it) tau_star of u0
    Au0_b       = Ax(u0) - b;
    u0_tau_star = norm( Au0_b(:) ) / (sqrt(d)*n_sigma);
    % print preliminaries of model/algorithm/initial iterate in the command window
    fprintf('\n\n\nIMAGE %s by TV-L2-DC MODEL solved by ADMM:\n',INV_PROB_NAME_LONG);
    if ( TV_L2_DC_taus_n == 1 ), itrs_debugs = 1; else, itrs_debugs = 0; 
        fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);    
    end
    % cycle over all the selected different tau values
    TV_L2_DC_ADMM_BEST_ISNR  = -10^5; 
    TV_L2_DC_ADMM_BEST_ISSIM = -10^5;
    for tau_i = 1:TV_L2_DC_taus_n
        tau = TV_L2_DC_taus(tau_i);
        % for the current value of the discrepancy parameter tau, 
        % compute the reconstructed image u_star(tau) by solving the model by ADMM        
        t0 = cputime;
        
        if ( TEST_REST_OR_INPT == 0 )
            [u_star,itrs] = REST_TV_L2_DC_ADMM(b,blur_k,n_sigma,tau,TV_L2_DC_ADMM_beta_t,TV_L2_DC_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,u_true);
        else
            [u_star,itrs] = INPT_TV_L2_DC_ADMM(b,M,n_sigma,tau,TV_L2_DC_ADMM_beta_t,TV_L2_DC_ADMM_beta_r,u0,itrs_th,itrs_rel_chg_th,itrs_debugs,u_true);
        end
        
        % compute/store the current reconstruction scalar performance measures
        TV_L2_DC_ADMM_CPU_TIMES(tau_i) = cputime - t0;
        TV_L2_DC_ADMM_ISNRS(tau_i)     = compute_snr(u_star,u_true) - u0_SNR;
        TV_L2_DC_ADMM_ISSIMS(tau_i)    = ssim(u_star,u_true) - u0_SSIM;
        TV_L2_DC_ADMM_ITRS(tau_i)      = itrs;
        Au_star_b                      = Ax(u_star) - b;
        TV_L2_DC_ADMM_TAU_STARS(tau_i) = norm( Au_star_b(:) ) / (sqrt(d)*n_sigma);
        % update the best reconstruction results in terms of accuracy measures ISNR and SSIM 
        if ( TV_L2_DC_ADMM_ISNRS(tau_i) > TV_L2_DC_ADMM_BEST_ISNR )
            TV_L2_DC_ADMM_BEST_ISNR         = TV_L2_DC_ADMM_ISNRS(tau_i);
            TV_L2_DC_ADMM_BEST_ISNR_u_star  = u_star;
            TV_L2_DC_ADMM_BEST_ISNR_tau_i   = tau_i;
        end
        if ( TV_L2_DC_ADMM_ISSIMS(tau_i) > TV_L2_DC_ADMM_BEST_ISSIM )
            TV_L2_DC_ADMM_BEST_ISSIM        = TV_L2_DC_ADMM_ISSIMS(tau_i);
            TV_L2_DC_ADMM_BEST_ISSIM_u_star = u_star;
            TV_L2_DC_ADMM_BEST_ISSIM_tau_i  = tau_i;
        end
        % print in the Matlab command window the current reconstruction results
        if ( TV_L2_DC_taus_n == 1 )
            fprintf('\n                              u0:  tau* = %6.4f,  [ SNR, SSIM] = [%6.3f , %6.4f]',u0_tau_star,u0_SNR,u0_SSIM);
        end
        fprintf('\ntau(%02d/%02d) = %6.4f:  ITS = %04d,  tau* = %6.4f,  [ISNR,ISSIM] = [%6.3f , %6.4f]',...
            tau_i,TV_L2_DC_taus_n,tau,itrs,TV_L2_DC_ADMM_TAU_STARS(tau_i),TV_L2_DC_ADMM_ISNRS(tau_i),TV_L2_DC_ADMM_ISSIMS(tau_i));
        % show in a Matlab figure the current reconstruction results
        if ( SHOW_ALL_RESTORED_IMAGES == 1 )
            if ( tau_i == 1 )
                figure('Name','Figure 81:  TV-L2-DC-ADMM:   tau-DEPENDENT RESULTS','NumberTitle','off');
                set(gcf,'Position',get(0,'ScreenSize'));
            end
            subplot(2,2,1)
            imshow(uint8(IMS_VIS_SF*u_true));
            title(sprintf('IM%02d (%d x %d):   u_{true}',im_id,h,w));
            subplot(2,2,3)
            if ( TEST_REST_OR_INPT == 0 )
                imshow(uint8(IMS_VIS_SF*b));
            else
                imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
            end
            title(sprintf('b = N(A u_{true}):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM));
            subplot(2,2,2)
            imshow(uint8(IMS_VIS_SF*u_star));
            title(sprintf('u^*(tau = %6.4f):  [ISNR,ISSIM] = [%6.3f,%6.4f]',tau,TV_L2_DC_ADMM_ISNRS(tau_i),TV_L2_DC_ADMM_ISSIMS(tau_i)));
            subplot(2,2,4)
            imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star) )); 
            title(sprintf('| u^*(tau = %6.4f) - u_{true} |',tau));
            pause(0.1);
        end
        % save on the hard disk the current reconstruction results
        if ( SAVE_ALL_RESTORED_IMAGES == 1 )
            imwrite(uint8(IMS_VIS_SF*u_star),sprintf('%s/%s_IM%02d_D8_TV_L2_DC_ADMM_A_REST_tau_%06.4f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,tau));
            imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-u_star)),sprintf('%s/%s_IM%02d_D8_TV_L2_DC_ADMM_B_AERR_tau_%06.4f.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,tau));
        end
    end
    % print in the Matlab command window the best reconstruction results in terms of ISNR and ISSIM
    TV_L2_DC_ADMM_BEST_ISNR_tau       = TV_L2_DC_taus(TV_L2_DC_ADMM_BEST_ISNR_tau_i);
    TV_L2_DC_ADMM_BEST_ISSIM_tau      = TV_L2_DC_taus(TV_L2_DC_ADMM_BEST_ISSIM_tau_i);
    TV_L2_DC_ADMM_BEST_ISNR_itrs      = TV_L2_DC_ADMM_ITRS(TV_L2_DC_ADMM_BEST_ISNR_tau_i);
    TV_L2_DC_ADMM_BEST_ISSIM_itrs     = TV_L2_DC_ADMM_ITRS(TV_L2_DC_ADMM_BEST_ISSIM_tau_i);
    TV_L2_DC_ADMM_BEST_ISNR_tau_star  = TV_L2_DC_ADMM_TAU_STARS(TV_L2_DC_ADMM_BEST_ISNR_tau_i);
    TV_L2_DC_ADMM_BEST_ISSIM_tau_star = TV_L2_DC_ADMM_TAU_STARS(TV_L2_DC_ADMM_BEST_ISSIM_tau_i);        
    fprintf('\nBEST ISNR and ISSIM RESULTS:');
    fprintf('\ntau(%02d/%02d) = %6.4f:  ITS = %04d,  tau* = %6.4f,   BEST ISNR   =  %6.3f',...
        TV_L2_DC_ADMM_BEST_ISNR_tau_i,TV_L2_DC_taus_n,TV_L2_DC_ADMM_BEST_ISNR_tau,...
        TV_L2_DC_ADMM_BEST_ISNR_itrs,TV_L2_DC_ADMM_BEST_ISNR_tau_star,TV_L2_DC_ADMM_BEST_ISNR);
    fprintf('\ntau(%02d/%02d) = %6.4f:  ITS = %04d,  tau* = %6.4f,   BEST ISSIM  =   %6.4f',...
        TV_L2_DC_ADMM_BEST_ISSIM_tau_i,TV_L2_DC_taus_n,TV_L2_DC_ADMM_BEST_ISSIM_tau,...
        TV_L2_DC_ADMM_BEST_ISSIM_itrs,TV_L2_DC_ADMM_BEST_ISSIM_tau_star,TV_L2_DC_ADMM_BEST_ISSIM);
    % show in a Matlab figure the best reconstruction results in terms of ISNR and ISSIM
    if ( SHOW_BEST_RESTORED_IMAGES == 1 )
        figure('Name','Figure 82:  TV-L2-DC-ADMM:  BEST RESULTS in terms of ISNR and ISSIM','NumberTitle','off');
        set(gcf,'Position',get(0,'ScreenSize'));
        subplot(2,3,1)
        imshow(uint8(IMS_VIS_SF*u_true));
        title(sprintf('IM%02d (%d x %d):   u_{true}',im_id,h,w));
        subplot(2,3,4)
        if ( TEST_REST_OR_INPT == 0 )
            imshow(uint8(IMS_VIS_SF*b));
        else
            imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
        end
        title(sprintf('b = N(A u_{true}):  [SNR,SSIM] = [%6.3f,%6.4f]',u0_SNR,u0_SSIM));
        subplot(2,3,2)
        imshow(uint8(IMS_VIS_SF*TV_L2_DC_ADMM_BEST_ISNR_u_star));
        title(sprintf('u^*(tau = %6.4f):  BEST ISNR = %6.3f',TV_L2_DC_ADMM_BEST_ISNR_tau,TV_L2_DC_ADMM_BEST_ISNR));
        subplot(2,3,5)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L2_DC_ADMM_BEST_ISNR_u_star) )); 
        title(sprintf('| u^*(tau = %6.4f) - u_{true} |',TV_L2_DC_ADMM_BEST_ISNR_tau));
        subplot(2,3,3)
        imshow(uint8(IMS_VIS_SF*TV_L2_DC_ADMM_BEST_ISSIM_u_star));
        title(sprintf('u^*(tau = %6.4f):  BEST ISSIM = %6.4f',TV_L2_DC_ADMM_BEST_ISSIM_tau,TV_L2_DC_ADMM_BEST_ISSIM));
        subplot(2,3,6)
        imshow(uint8( AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L2_DC_ADMM_BEST_ISSIM_u_star) )); 
        title(sprintf('| u^*(tau = %6.4f) - u_{true} |',TV_L2_DC_ADMM_BEST_ISSIM_tau));
    end
    % save on the hard disk the best reconstruction results in terms of ISNR and ISSIM
    if ( SAVE_BEST_RESTORED_IMAGES == 1 )
        imwrite(uint8(IMS_VIS_SF*TV_L2_DC_ADMM_BEST_ISNR_u_star),sprintf('%s/%s_IM%02d_D8_TV_L2_DC_ADMM_A_REST_tau_%06.4f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L2_DC_ADMM_BEST_ISNR_tau));
        imwrite(uint8(IMS_VIS_SF*TV_L2_DC_ADMM_BEST_ISSIM_u_star),sprintf('%s/%s_IM%02d_D8_TV_L2_DC_ADMM_A_REST_tau_%06.4f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L2_DC_ADMM_BEST_ISSIM_tau));
        imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L2_DC_ADMM_BEST_ISNR_u_star)),sprintf('%s/%s_IM%02d_D8_TV_L2_DC_ADMM_B_AERR_tau_%06.4f_BEST_ISNR.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L2_DC_ADMM_BEST_ISNR_tau));
        imwrite(uint8(AEI_VIS_SF*IMS_VIS_SF*abs(u_true-TV_L2_DC_ADMM_BEST_ISSIM_u_star)),sprintf('%s/%s_IM%02d_D8_TV_L2_DC_ADMM_B_AERR_tau_%06.4f_BEST_ISSIM.png',SAVE_FOLDER,INV_PROB_NAME_SHORT,im_id,TV_L2_DC_ADMM_BEST_ISSIM_tau));
    end
    % show in a Matlab figure the performance graphs
    % show...only if we have selected at least two different tau values
    taus_min = min(TV_L2_DC_taus); taus_max = max(TV_L2_DC_taus);
    if ( taus_max > taus_min )
        if ( SHOW_PERFORMANCE_GRAPHS_VS_MU == 1 )
            figure('Name','Figure 83:  TV-L2-DC-ADMM:  PERFORMANCE GRAPHS  versus  TAU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TV-L2-DC-ADMM:    ACCURACY MEASURES  versus  \tau');
            yyaxis left;
            plot_y_min = min(TV_L2_DC_ADMM_ISNRS); plot_y_max = max(TV_L2_DC_ADMM_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_DC_taus,TV_L2_DC_ADMM_ISNRS,'-o','markersize',3);
            xlabel('\tau'); ylabel('ISNR(\tau)  [dB]'); 
            axis([taus_min,taus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L2_DC_ADMM_ISSIMS); plot_y_max = max(TV_L2_DC_ADMM_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_DC_taus,TV_L2_DC_ADMM_ISSIMS,'-o','markersize',3);
            xlabel('\tau'); ylabel('ISSIM(\tau)');
            axis([taus_min,taus_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TV-L2-DC-ADMM:    EFFICIENCY MEASURES  versus  \tau');
            yyaxis left;
            plot_y_min = min(TV_L2_DC_ADMM_ITRS); plot_y_max = max(TV_L2_DC_ADMM_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_DC_taus,TV_L2_DC_ADMM_ITRS,'-o','markersize',3);
            xlabel('\tau'); ylabel('ITERS(\tau)');
            axis([taus_min,taus_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L2_DC_ADMM_CPU_TIMES); plot_y_max = max(TV_L2_DC_ADMM_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_DC_taus,TV_L2_DC_ADMM_CPU_TIMES,'-o','markersize',3);
            xlabel('\tau'); ylabel('CPU-TIME(\tau)  [secs]');
            axis([taus_min,taus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_TAU_STAR_VS_MU == 1 )
            figure('Name','Figure 84:  TV-L2-DC-ADMM:  TAU*  versus  TAU','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));            
            plot_y_min = min(TV_L2_DC_ADMM_TAU_STARS); plot_y_max = max(TV_L2_DC_ADMM_TAU_STARS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_DC_taus,TV_L2_DC_ADMM_TAU_STARS,'-om','markersize',3);
            xlabel('\tau'); ylabel('\tau^*(\tau)');
            title('TV-L2-DC-ADMM:    TAU*  versus  \tau');
            axis([taus_min,taus_max,plot_y_min,plot_y_max]);
        end
        if ( SHOW_PERFORMANCE_GRAPHS_VS_TAU_STAR == 1 )
            tau_stars_min = min(TV_L2_DC_ADMM_TAU_STARS); tau_stars_max = max(TV_L2_DC_ADMM_TAU_STARS);
            figure('Name','Figure 85:  TV-L2-DC-ADMM:  PERFORMANCE GRAPHS  versus  TAU*(TAU)','NumberTitle','off');
            set(gcf,'Position',get(0,'ScreenSize'));
            subplot(2,1,1)
            title('TV-L2-DC-ADMM:    ACCURACY MEASURES  versus  \tau^*(\tau)');
            yyaxis left;
            plot_y_min = min(TV_L2_DC_ADMM_ISNRS); plot_y_max = max(TV_L2_DC_ADMM_ISNRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_DC_ADMM_TAU_STARS,TV_L2_DC_ADMM_ISNRS,'-o','markersize',3);
            xlabel('\tau^*(\tau)'); ylabel('ISNR(\tau^*(\tau))  [dB]'); 
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L2_DC_ADMM_ISSIMS); plot_y_max = max(TV_L2_DC_ADMM_ISSIMS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_DC_ADMM_TAU_STARS,TV_L2_DC_ADMM_ISSIMS,'-o','markersize',3);
            xlabel('\tau^*(\tau)'); ylabel('ISSIM(\tau^*(\tau))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            subplot(2,1,2)
            title('TV-L2-DC-ADMM:    EFFICIENCY MEASURES  versus  \tau^*(\tau)');
            yyaxis left;
            plot_y_min = min(TV_L2_DC_ADMM_ITRS); plot_y_max = max(TV_L2_DC_ADMM_ITRS);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot_y_max = plot_y_max + 0.05 * (plot_y_max-plot_y_min);
            plot(TV_L2_DC_ADMM_TAU_STARS,TV_L2_DC_ADMM_ITRS,'-o','markersize',3);
            xlabel('\tau^*(\tau)'); ylabel('ITERS(\tau^*(\tau))');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
            yyaxis right;
            plot_y_min = min(TV_L2_DC_ADMM_CPU_TIMES); plot_y_max = max(TV_L2_DC_ADMM_CPU_TIMES);
            plot_y_rng = max(plot_y_max - plot_y_min,0.001);
            plot_y_min = plot_y_min - 0.04 * plot_y_rng; plot_y_max = plot_y_max + 0.04 * plot_y_rng;
            plot(TV_L2_DC_ADMM_TAU_STARS,TV_L2_DC_ADMM_CPU_TIMES,'-o','markersize',3);
            xlabel('\tau^*(\tau)'); ylabel('CPU-TIME(\tau^*(\tau))  [secs]');
            axis([tau_stars_min,tau_stars_max,plot_y_min,plot_y_max]);
        end
    end
end

%--------------------------------------------------------------------------


% finally, stop (eventually, if you chose to start it) CPU time profiling
if ( PROFILE_CPU_TIMES == 1 )
    profile viewer;  % view the CPU time profiler results
    profile off;     % stop CPU time profiler
end

