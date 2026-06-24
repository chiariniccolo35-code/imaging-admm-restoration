function [u] = REST_TIK_L2_U_DIR(b,blur_k,mu)
%
% Compute approximate solutions (that is, minimizers) of the
% unconstrained TIK-L2 variational model: 
%
% u*(mu) = argmin J(u;mu),   J(u;mu) = (1/2)||D u||_2^2 + (mu/2)||K u - b||_2^2  
%             u
%
% for the non-blind restoration, that is denoising + deblurring (non-blind 
% means that the blur operator is assumed to be known), of grayscale images 
% using a "direct" solver by 2D FFT.
%
% We can compute the solution in closed-form by solving the linear system:
%
% (DhT Dh + DvT Dv+ mu KT K)x = mu KT b ... see the course slides.
%
% Since we are assuming space-invariant blur corruptions, so that the 
% blur operator is a 2D convolution operator, we can solve the linear 
% system by fast 2D transforms, according to chosen boundary conditions.
% We assume periodic boundary conditions -> direct solution by 2D FFT.

% INPUTS:
    % inputs defining the unconstrained TIK-L2 cost function to be minimized:
        % b         --> observed corrupted (blurred and noisy) image
        % blur_k    --> blur PSF (kernel of spatial 2D convolution)
        % mu        --> regularization parameter (positive real scalar)
        
% We remark the following:
%   - all images are stored as matrices, that is they are not vectorized!
%   - periodic boundary conditions are assumed for the unknown clean image 
%     to restore, that is both the blur and the first-order derivatives 
%     near the image boundaries are computed according to this assumption!
%     This assumption allows to diagonalize the blur matrix K and the 
%     first-order horizontal and vertical partial derivative matrices
%     Dh and Dv by means of the 2D Discrete Fourier Transform (DFT)
    
% extract image dimensions
[h,w] = size(b); % height and width, in pixels

% Fourier-diagonalize the matrices Dh, Dv, K
Dh_DFT  = psf2otf([1,-1],[h,w]);
Dv_DFT  = psf2otf([1;-1],[h,w]);
K_DFT   = psf2otf(blur_k,[h,w]);

% compute the Fourier-diagonalized coefficient matrix
CM_DFT  =      conj(Dh_DFT) .* Dh_DFT + ...
               conj(Dv_DFT) .* Dv_DFT + ...
          mu * conj(K_DFT)  .* K_DFT;

% compute the Fourier-transformed right-hand side
rhs_DFT = mu * conj(K_DFT) .* fft2(b);

% compute the linear system unique solution
u       = real(ifft2( rhs_DFT ./ CM_DFT ));


end


