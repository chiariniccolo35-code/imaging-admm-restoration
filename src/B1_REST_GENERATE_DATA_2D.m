close all; clear all; clc;


SHOW_CORRUPTIONS      = 1;
SHOW_NOISE_HISTOGRAMS = 1;


% -------------------------------------------------------------------------
% SELECT THE ORIGINAL IMAGE TO BE TESTED
im_id = 5; % index of the original image to be tested (see the switch-case below)
switch im_id    
    % PURELY PIECEWISE CONSTANT IMAGES
    case 0 % synthetic binary square
        im_file     = './test_images/00_square.png';        
        u_true      = double(imread(im_file))/255;
    case 1 % synthetic binary rectangles
        im_file     = './test_images/01_rectangles.png';        
        u_true      = double(imread(im_file))/255;
    case 2 % checkboard
        im_file     = './test_images/02_checkboard.png';        
        u_true      = double(imread(im_file))/255;
    case 3 % checkboard fine
        im_file     = './test_images/03_checkboard_finer.png';        
        u_true      = double(imread(im_file))/255;
    case 4 % qrcode
        im_file     = './test_images/04_qrcode.png';
        u_true      = double(imread(im_file))/255;
    case 5 % qrcode_finer
        im_file     = './test_images/05_qrcode_finer.png';
        u_true      = double(imread(im_file))/255;
    case 6 % geometric_A
        im_file     = './test_images/06_geometric_A.png';
        u_true      = double(imread(im_file))/255;
    case 7 % geometric_B
        im_file     = './test_images/07_geometric_B.png';
        u_true      = double(imread(im_file))/255;        
    case 8 % head
        im_file     = './test_images/08_head.png';
        u_true      = double(imread(im_file))/255;
    % PREVALENTLY PIECEWISE CONSTANT IMAGES
    case 10 % satellite
        im_file     = './test_images/10_satellite.png';
        u_true      = double(imread(im_file))/255;
    case 11 % cameraman
        im_file     = './test_images/11_cameraman.png';
        u_true      = double(imread(im_file))/255;
    case 12 % brain_section
        im_file     = './test_images/12_brain_section.png';
        u_true      = double(imread(im_file))/255;
    case 13 % mri
        im_file     = './test_images/13_mri.png';
        u_true      = double(imread(im_file))/255;
    % PIECEWISE SMOOTH IMAGES
    case 21 % peppers
        im_file     = './test_images/21_peppers256.png';
        u_true      = double(imread(im_file))/255;
    case 22 % Lena
        im_file     = './test_images/22_lena256.png';
        u_true      = double(imread(im_file))/255;
    case 23 % Elaine
        im_file     = './test_images/23_elaine.png';
        u_true      = double(imread(im_file))/255;
    case 24 % girlface
        im_file     = './test_images/24_girlface.png';
        u_true      = double(imread(im_file))/255;
    case 25 % butterfly
        im_file     = './test_images/25_butterfly.png';
        u_true      = double(imread(im_file))/255;
    % PURELY SMOOTH IMAGES
    case 30 % synthetic 2D synusoids
        im_file     = './test_images/30_sinusoids.png';
        u_true      = double(imread(im_file))/255;
    % SOME "NATURAL" IMAGES
    case 40 % bridge 
        im_file     = './test_images/40_bridge256.png';
        u_true      = double(imread(im_file))/255;
    case 41 % boats 
        im_file     = './test_images/41_boats256.png';
        u_true      = double(imread(im_file))/255;
    case 42 % hill 
        im_file     = './test_images/42_hill256.png';
        u_true      = double(imread(im_file))/255;
end

% extract image dimensions
[h,w] = size(u_true); % height and width, in pixels
d     = h * w;        % total number of pixels

% -------------------------------------------------------------------------
% SET THE BLUR PARAMETERS, THEN GENERATE AND STORE THE BLUR KERNEL 
% (or PSF) blur_k, FINALLY BLUR THE ORIGINAL IMAGE u_true --> Au_true

    % set the blur parameters then generate and store the blur kernel, blur_k
    blur_type_id = 3; % index of the blur type (see the switch-case below)
    switch blur_type_id
        case 0 % no blur
            blur_k             = 1; % kernel
            blur_type_descr    = 'NO'; % type description            
        case 1 % average
            blur_r             = 1; % radius (pixels)
            blur_k             = fspecial('average',(1 + 2 * blur_r)); % kernel 
            blur_bc_type       = 0; % boundary conditions type (0->periodic)
            blur_type_descr    = 'AVER'; % type description          
        case 2 % disk
            blur_r             = 3; % radius (pixels)
            blur_k             = fspecial('disk',blur_r); % kernel
            blur_bc_type       = 0; % boundary conditions type (0->periodic)
            blur_type_descr    = 'DISK'; % type description
        case 3 % Gaussian
            blur_r             = 2; % radius (pixels) ...the band is 2 * blur_r + 1
            blur_s             = 1; % standard deviation
            blur_k             = fspecial('gaussian',(1 + 2 * blur_r),blur_s); % kernel
            blur_bc_type       = 0; % boundary conditions type (0->periodic)
            blur_type_descr    = 'GAUS'; % type description
        case 4 % motion
            blur_l             = 13; % length (pixels) 
            blur_a             = 30; % angle
            blur_k             = fspecial('motion',blur_l,blur_a); % kernel
            blur_bc_type       = 0; % boundary conditions type (0->periodic)
            blur_type_descr    = 'MOTION'; % type description
    end
    % normalize the blur kernel such that it integrates to 1 
    % (actually, it is done by fspecial!)
    blur_k = blur_k / sum(blur_k(:)); 
    % starting from the blur kernel (or Point Spread Function, PSF) blur_k, 
    % compute the so-called blur Optical Transfer Function (OTF), which 
    % mathematically represents the Discrete Fourier Transform (DFT) 
    % diagonalization of the blur convolution matrix K, hence denoted K_DFT, 
    % then compute (in Fourier domain) and store the blurred image Au_true
    if ( blur_type_id == 0 ) % no blur, K is the identity matrix
        K_DFT   = psf2otf(blur_k,size(u_true));
        Au_true = u_true;
    else % blur, K is a blur convolution matrix
        switch blur_bc_type
            case 0 % periodic -> Fourier
                K_DFT   = psf2otf(blur_k,size(u_true));
                Au_true = real( ifft2( K_DFT .* fft2(u_true) ) );
            case 1 % Neumann homogeneous
                % ...              
        end
    end

% -------------------------------------------------------------------------
% SET NOISE PARAMETERS, THEN "ADD" NOISE TO THE BLURRED IMAGE Au_true --> b
noise_type_id = 2; % index of the noise type (see the switch-case below)
switch noise_type_id
    case 0 % no noise
        n_mean          = 0/255;   % mean (intended for images in [0,1])
        n_sigma         = 0/255;   % stdv (intended for images in [0,1])
        n_type_descr    = 'NO';    % type description
        n_realiz        = zeros(h,w);
        b               = Au_true;
    case 1 % additive white Gaussian noise (AWGN)
        n_mean          = 0/255;   % mean (intended for images in [0,1])
        n_sigma         = 15/255;  % stdv (intended for images in [0,1])
        n_type_descr    = 'AWG';   % type description
        n_realiz        = n_mean + n_sigma * randn(h,w);
        b               = Au_true + n_realiz;
    case 2 % additive white uniform noise (AWUN) 
           % (given n_mean and n_sigma -> uniform pdf in [n_mean-sqrt(3)n_sigma,n_mean+sqrt(3)n_sigma])
        n_mean          = 0/255;   % mean (intended for images in [0,1])
        n_sigma         = 30/255;  % stdv (intended for images in [0,1])
        n_type_descr    = 'AWU';   % type description
        n_realiz        = n_mean + (sqrt(3) * n_sigma) * ( 2 * (rand(h,w) - 0.5) );
        b               = Au_true + n_realiz;
    case 3 % additive white Laplacian noise (AWLN)
        n_mean          = 0/255;   % mean (intended for images in [0,1])
        n_sigma         = 30/255;  % stdv (intended for images in [0,1])
        n_type_descr    = 'AWL';   % type description
        n_sp            = n_sigma / sqrt(2); % scale parameter
        n_U             = rand(h,w) - 0.5;
        n_realiz        = n_mean - n_sp * ( sign(n_U) .* log(1 - 2 * abs(n_U)) );
        b               = Au_true + n_realiz;
    case 4 % additive white generalized Gaussian noise (AWGGN)
        n_mean          = 0/255;  % mean (intended for images in [0,1])
        n_sigma         = 30/255; % stdv (intended for images in [0,1])
        n_q             = 0.3;    % shape parameter 
        n_type_descr    = 'AWGG'; % type description
        n_realiz        = generate_WGG_realization(h,w,n_mean,n_sigma,n_q);
        b               = Au_true + n_realiz;
    case 5 % multiplicative white Gaussian noise (MWGN)
        n_mean          = 1;        % mean
        n_sigma         = 20/127.5; % stdv
        n_type_descr    = 'MWG';    % type description
        n_realiz        = n_mean + n_sigma * randn(h,w);
        b               = Au_true .* n_realiz;
    case 6 % multiplicative white uniform noise (MWUN)
        n_mean          = 1;        % mean
        n_sigma         = 20/127.5; % stdv
        n_type_descr    = 'MWU';    % type description
        n_realiz        = n_mean + (sqrt(3) * n_sigma) * ( 2 * (rand(h,w) - 0.5) );
        b               = Au_true .* n_realiz;
    case 7 % multiplicative white Laplacian noise (MWLN)
        n_mean          = 1;        % mean
        n_sigma         = 20/127.5; % stdv
        n_type_descr    = 'MWL';    % type description
        n_sp            = n_sigma / sqrt(2); % scale parameter
        n_U             = rand(h,w) - 0.5;
        n_realiz        = n_mean - n_sp * ( sign(n_U) .* log(1 - 2 * abs(n_U)) );
        b               = Au_true .* n_realiz;
    case 8 % multiplicative white generalized Gaussian noise (MWGGN)
        n_mean          = 1;        % mean
        n_sigma         = 20/127.5; % stdv
        n_q             = 0.3;    % shape parameter 
        n_type_descr    = 'MWGG';   % type description
        n_realiz        = generate_WGG_realization(h,w,n_mean,n_sigma,n_q);
        b               = Au_true .* n_realiz;
    case 9 % impulsive: salt & pepper
        n_p             = 0.4; % probability  for a pixel to be noise-corrupted
        n_type_descr    = 'ISP';
        n_ph            = n_p/2;
        n_P             = rand(h,w);
        b               = Au_true;
        b(n_P<=n_ph)    = 0;
        b(n_P>(1-n_ph)) = 1;
        n_sigma         = std(b(:) - Au_true(:));
    case 10 % impulsive: random-valued
        n_p             = 0.1; % probability  for a pixel to be noise-corrupted
        n_type_descr    = 'IRV';
        n_P             = rand(h,w);
        n_P1            = rand(h,w);
        b               = Au_true;
        b(n_P<=n_p)     = n_P1(n_P<=n_p);
        n_sigma         = std(b(:) - Au_true(:));
    case 11 % Poisson
        POIS_fact       = 20;
        u_true          = POIS_fact * u_true;
        Au_true         = POIS_fact * Au_true;
        n_type_descr    = 'POIS';
        b               = poissrnd(Au_true);
        n_sigma         = std(b(:) - Au_true(:));
end

% depending on the selected noise type, set the scaling
% factor for visualization/saving of all the images
if ( noise_type_id <= 10 )
    IMS_VIS_SF = 255;
elseif ( noise_type_id == 11 )
    IMS_VIS_SF = 255 / (1.2*POIS_fact);
end

if SHOW_CORRUPTIONS
    figure(1)
    set(gcf,'Position',get(0,'ScreenSize'));
    subplot(2,3,[1,4])
    imshow(uint8(IMS_VIS_SF*u_true));
    title(sprintf('ORIGINAL (%d x %d)',h,w));
    subplot(2,3,2)
    imshow(uint8(IMS_VIS_SF*Au_true));
    title(sprintf('CORRUPTED by %s BLUR',blur_type_descr));
    subplot(2,3,3)
    imshow(uint8(IMS_VIS_SF*b));
    title(sprintf('CORRUPTED by %s BLUR and %s NOISE',blur_type_descr,n_type_descr));
    subplot(2,3,5)
    [blur_k_h,blur_k_w] = size(blur_k);
    enl     = 1;
    blur_k_2   = zeros(blur_k_h+2*enl,blur_k_w+2*enl);
    blur_k_2((1+enl):(end-enl),(1+enl):(end-enl)) = blur_k;
    image(255*blur_k_2/max(blur_k_2(:)));
    axis equal;
    axis tight;
    title(sprintf('BLUR KERNEL (PSF): %s %dx%d',blur_type_descr,size(blur_k,1),size(blur_k,2)));
    subplot(2,3,6)
    imshow(uint8(127.5+IMS_VIS_SF*(b-Au_true)));
    title(sprintf('%s NOISE',n_type_descr));
end

if SHOW_NOISE_HISTOGRAMS
    if ( ( noise_type_id >=1) && ( noise_type_id <= 8) )
        % compute histogram of noise realization
        n_H_edges_n       = 100;
        n_H_edges_min     = n_mean - 4 * n_sigma;
        n_H_edges_max     = n_mean + 4 * n_sigma;
        n_H_edges         = linspace(n_H_edges_min,n_H_edges_max,n_H_edges_n);
        n_H_centers       = 0.5 * ( n_H_edges(1:(n_H_edges_n-1)) + n_H_edges(2:n_H_edges_n) );
        n_H_bins_n        = n_H_edges_n - 1;
        n_H_bins_size     = n_H_edges(2) - n_H_edges(1);        
        n_H               = histcounts( n_realiz , n_H_edges );
        n_H               = n_H(1:n_H_bins_n) / sum( n_H(1:n_H_bins_n) );
        n_H               = n_H / n_H_bins_size;
        % compute theoretical noise pdf
        n_pdf_xs = linspace(n_H_edges_min,n_H_edges_max,10000);
        switch noise_type_id 
            case {1,5} % gaussian
                n_pdf = gauss_pdf_1D(n_mean,n_sigma,n_pdf_xs);
            case {2,6} % uniform
                n_pdf = uniform_pdf_1D(n_mean,n_sigma,n_pdf_xs);
            case {3,7} % laplacian
                n_pdf = laplace_pdf_1D(n_mean,n_sigma,n_pdf_xs);
            case {4,8} % generalized Gaussian
                n_pdf = generalized_gauss_pdf_1D(n_mean,n_sigma,n_q,n_pdf_xs);
        end
        % create figure
        figure(2)
        hold on;
        bar(n_H_centers,n_H,'b');
        plot(n_pdf_xs,n_pdf,'r','linewidth',2);
        axis([n_H_edges_min,n_H_edges_max,0,1.05*max([max(n_H),max(n_pdf)])]);
        legend('histogram','theoretical pdf');
        xlabel('x');
        ylabel('p(x)');
        title(sprintf('HISTOGRAM and THEORETICAL PDF of %s NOISE REALIZATION',n_type_descr));        
    end
    
end



