function [u] = INPT_TIK_L2_U_DIR(b,M,mu)
%
% Compute approximate solutions (that is, minimizers) of the
% unconstrained TIK-L2 variational model: 
%
% u*(mu) = argmin J(u;mu),   J(u;mu) = (1/2)||D u||_2^2 + (mu/2)||S u - b||_2^2  
%             u
%
% for the non-blind inpainting, that is denoising + holes filling (non-blind 
% means that the selection/masking operator is assumed to be known)
%
% We can compute the solution in closed-form by solving the linear system:
%
% (DhT Dh + DvT Dv+ mu S)x = mu b ... see the course slides.
%
% We create (in sparse format) the coefficient matrix and the right-hand 
% side vector, then we solve the linear system directly by backslash \ ,
% which will use a (sparse) Cholesky solver, as the coefficient matrix 
% of the linear system above is symmetric positive definite and sparse

% INPUTS:
    % inputs defining the unconstrained TIK-L2 cost function to be minimized:
        % b         --> observed corrupted (masked and noisy) image
        % M         --> inpainting (binary) mask image
        % mu        --> regularization parameter (positive real scalar)
        
% We remark the following:
%   - all images are stored as matrices, that is they are not vectorized!
%   - periodic boundary conditions are assumed for the unknown clean image 
%     to restore, that is the first-order derivatives near the image 
%     boundaries are computed according to this assumption!

% extract image dimensions
[h,w] = size(b); % height and width, in pixels
d     = h * w;   % total number of pixels

% generate (in sparse format) the matrices Dh, Dv and S
    % Dh
    dh           = ones(w,1);
    Dh_1d        = spdiags([-dh,dh],[0,1],w,w);
    Dh_1d(end,1) = 1;
    Dh           = kron( Dh_1d , speye(h) );
    % Dv
    dv           = ones(h,1);
    Dv_1d        = spdiags([-dv,dv],[0,1],h,h);
    Dv_1d(end,1) = 1;
    Dv           = kron( speye(w) , Dv_1d );
    % S
    S            = spdiags(M(:),0,d,d);

% compute (in sparse format) the coefficient matrix
CM = Dh'*Dh + Dv'*Dv + mu * S;

% compute (in sparse format) the right-hand side
rhs = mu * sparse(b(:));

% compute the linear system unique solution
u  = reshape( full(CM \ rhs) , h,w);


end


