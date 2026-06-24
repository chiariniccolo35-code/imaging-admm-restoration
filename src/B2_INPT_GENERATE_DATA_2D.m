%close all; clear all; clc;


SHOW_CORRUPTIONS      = 1;
SHOW_NOISE_HISTOGRAMS = 0;


% -------------------------------------------------------------------------
% SELECT THE ORIGINAL IMAGE TO BE TESTED
im_id = 43; % index of the original image to be tested (see the switch-case below)
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

    case 43 % Peppers
        im_file     = './test_images/peppers_true.png';
        u_true      = double(imread(im_file)) / 255;

        im_file_1   = './test_images/peppers_degr.png';
        u_degr      = double(imread(im_file_1)) / 255;

    case 44 % QR code
        im_file     = './test_images/qr_code_true.png';
        u_true      = double(imread(im_file)) / 255;

        im_file_1   = './test_images/qr_code_degr.png';
        u_degr      = double(imread(im_file_1)) / 255;
        

    case 45 % Sinusoid
        im_file     = './test_images/sinusoid_true.png';
        u_true      = double(imread(im_file)) / 255;

        im_file_1     = './test_images/sinusoid_degr.png';
        u_degr      = double(imread(im_file_1)) / 255;

end

% extract image dimensions
[h,w] = size(u_true); % height and width, in pixels
d     = h * w;        % total number of pixels

% -------------------------------------------------------------------------
% SET THE (SYNTHETIC) INPAINTING MASK PARAMETERS, THEN GENERATE AND STORE  
% THE INPAINTING_MASK IMAGE, FINALLY MASK THE ORIGINAL IMAGE u_true --> Au_true

    % set the masking parameters then generate and store the mask image, M
    mask_type_id = 4; % index of the mask type (see the switch-case below)
    switch mask_type_id
        case 0 % no masking
            M                  = ones(h,w); % mask image
            mask_type_descr    = 'NO'; % type description            
        case 1 % random points
            P_mask             = 0.4; % probability  for a pixel to be masked
            M                  = rand(h,w);
            M(M<=P_mask)       = 0;
            M(M>P_mask)        = 1;
            mask_type_descr    = 'RPS';            
        case 2 % random horizontal strips
            strips_n           = 15; % number of strips
            strips_ht          = 1; % half-thickness of strips (in pixels)
            js                 = round(1 + (h-1) * rand(1,strips_n));
            M                  = ones(h,w);
            for strip_i = 1:strips_n
                j     = js(strip_i);
                j_min = max(j-strips_ht,1);
                j_max = min(j+strips_ht,h);
                M(j_min:j_max,:) = 0;
            end
            mask_type_descr    = 'RHSS'; % type description
        case 3 % random vertical strips
            strips_n           = 15; % number of strips
            strips_ht          = 1; % half-thickness of strips (in pixels)
            is                 = round(1 + (w-1) * rand(1,strips_n));
            M                  = ones(h,w);
            for strip_i = 1:strips_n
                i     = is(strip_i);
                i_min = max(i-strips_ht,1);
                i_max = min(i+strips_ht,w);
                M(:,i_min:i_max) = 0;
            end
            mask_type_descr    = 'RVSS'; % type description
        case 4 % super-resolution
            sr_factor          = 2;

            M                  = zeros(h, w);
            M( 1:sr_factor:h,1:sr_factor:h ) = 1;  

            mask_type_descr    = sprintf('SR%d',sr_factor); % type description
    end

    % mask the original image
      %
      Au_true = M .* u_true;

% -------------------------------------------------------------------------
% SET NOISE PARAMETERS, THEN "ADD" NOISE TO THE ORIGINAL IMAGE Au_true --> b
noise_type_id = 11; % index of the noise type (see the switch-case below)
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
        n_sigma         = 10/255;  % stdv (intended for images in [0,1])
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
        n_p             = 0.05; % probability  for a pixel to be noise-corrupted
        n_type_descr    = 'ISP';
        n_ph            = n_p/2;
        n_P             = rand(h,w);
        b               = Au_true;
        b(n_P<=n_ph)    = 0;
        b(n_P>(1-n_ph)) = 1;
        n_sigma         = std(b(M>0) - Au_true(M>0));
    case 10 % impulsive: random-valued
        n_p             = 0.1; % probability  for a pixel to be noise-corrupted
        n_type_descr    = 'IRV';
        n_P             = rand(h,w);
        n_P1            = rand(h,w);
        b               = Au_true;
        b(n_P<=n_p)     = n_P1(n_P<=n_p);
        n_sigma         = std(b(M>0) - Au_true(M>0));

    case 11 % LANZA EXERCISE
        M1                                = M;
        M1(1:sr_factor:h, 1:sr_factor:w)  = u_degr;
        b                                 = M1;
        n_sigma                           = std( b(M>0 ) - Au_true( M>0) );
        n_type_descr                      ='AWGG';
        
    case 12 % Poisson
        POIS_fact       = 20;
        u_true          = POIS_fact * u_true;
        Au_true         = POIS_fact * Au_true;
        n_type_descr    = 'POIS';
        b               = poissrnd(Au_true);
        n_sigma         = std(b(M>0) - Au_true(M>0));
     
end
% re-set to zero the inpainting mask pixels
%b(M==0) = 0;

% depending on the selected noise type, set the scaling
% factor for visualization/saving of all the images
if ( noise_type_id <= 11 )
    IMS_VIS_SF = 255;
elseif ( noise_type_id == 12 )
    IMS_VIS_SF = 255 / (1.2*POIS_fact);
end

% choose the color for visualization of inpainting mask
M_COL = [255;0;0];
    
if SHOW_CORRUPTIONS
    figure(1)
    set(gcf,'Position',get(0,'ScreenSize'));
    subplot(2,3,[1,4])
    imshow(uint8(IMS_VIS_SF*u_true));
    title(sprintf('ORIGINAL (%d x %d)',h,w),'Interpreter','Latex');
    subplot(2,3,2)
    imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*Au_true,M,M_COL) ));
    title(sprintf('CORRUPTED by %s MASK',mask_type_descr),'Interpreter','Latex');
    subplot(2,3,3)
    imshow(uint8( MASK_IMAGE_COL(IMS_VIS_SF*b,M,M_COL) ));
    title({
    sprintf('CORRUPTED by \\\\ %s MASK', mask_type_descr)
    sprintf('and %s NOISE', n_type_descr)
    }, 'interpreter', 'Latex');
    subplot(2,3,5)
    imshow(uint8(255*M));
    title(sprintf('%s BINARY INPAINTING MASK',mask_type_descr),'Interpreter','Latex');
    subplot(2,3,6)
    imshow(uint8( MASK_IMAGE_COL(127.5+IMS_VIS_SF*(b-Au_true),M,M_COL) ));
    title(sprintf('%s NOISE',n_type_descr),'interpreter','Latex');
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



