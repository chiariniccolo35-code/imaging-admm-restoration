clear all; close all; clc;

%in pratica per capire il noise abbiamo una sequenza di immagini prese di
%ste casette nella cartella SEQUENCE


% load from files the sequence of degraded images, 
% storing them in a pre-allocated 3d tensor  
ims_file_inds = 1:200;
ims_n         = numel(ims_file_inds);

% load the first image in order to read the images size
im           = imread(sprintf('./SEQUENCE/im_%03d.png',ims_file_inds(1)));
[h,w]        = size(im);
d            = h * w;

% pre-allocate the 3d tensor
ims_seq      = zeros(h,w,ims_n);

% load and store all the images
for im_i = 1:ims_n
    %im_file_ind       = ims_file_inds(im_i);
    im                = imread(sprintf('./SEQUENCE/im_%03d.png',im_i));
    ims_seq(:,:,im_i) = im;
    figure(1)
    imshow(im);
    title(sprintf('loading and storing image %d',im_i));
end

% compute, for each pixel position (i,j), the sample mean (mu_i,j) 
% and the sample standard deviation (sigma_i,j) of the sequence 
% of intensities (time series) that the pixel assume in the images

    % pre-allocate the mu and sigma matrices where
    % we will store all the mu_i,j and sigma_i,j values
    im_means = zeros(h,w);
    im_stdvs = zeros(h,w);

    % compute and store all the mu_i,j and sigma_i,j values, ciò che
    % andiamoa a calcolare è la media e la deviazione standard della
    % sequenza di 200 immagini pixel per pixel.
    for i = 1:h
        for j = 1:w
            time_series = ims_seq(i,j,:);
            im_means(i,j) = sum( time_series ) / ims_n;
            im_stdvs(i,j) = sqrt( sum( ( time_series - im_means(i,j) ).^2 ) / (ims_n - 1) ); % stimatore corretto della deviazione standard ha questa   
                                                                                             % espressione con n-1
        end
    end

% estimate the "best" sigma versus mu fitting curves, TEORIA SU NOTE
% according to the three noise models considered 
% (actually, we don't need to do that for the Poisson case):

    % AWGGN noise:  sigma = q (constant model)
      AWGG_q = sum( im_stdvs(:) ) / d;     % d = h * w;

    % MWGGN noise:  sigma = m mu (linear strictly increasing model) 
      MWGG  _m = sum( im_means(:) .* im_stdvs(:) ) / sum( im_means(:).^2 );

% compute the RMSR (Root Mean Square Residual) for the three models
AWGG_RMSR = sqrt( sum((im_stdvs(:) - AWGG_q).^2) / d );%residuo come scritto nelle note, sia uqi che sotto
MWGG_RMSR = sqrt( sum((im_stdvs(:) - MWGG_m * im_means(:)).^2) / d );
POIS_RMSR = sqrt( sum((im_stdvs(:) - sqrt(im_means(:))).^2) / d );%anche qui di base

% print in the command window the three RMSR 
% just computed for the three noise models
fprintf('\n');
fprintf('\nAWGG NOISE RMSR = %.3f',AWGG_RMSR);
fprintf('\nMWGG NOISE RMSR = %.3f',MWGG_RMSR);
fprintf('\nPOIS NOISE RMSR = %.3f',POIS_RMSR);
fprintf('\n');

% sort (in ascending order) the three RMSR ==> "best fitting" noise model
[RMSRS_sorted,RMSRS_sorted_inds] = sort([AWGG_RMSR,MWGG_RMSR,POIS_RMSR]);

% print in the command window the "best fitting" model (smallest RMSR)
switch RMSRS_sorted_inds(1)
    case 1 % AWGG
        fprintf('\nTHE BEST FITTING MODEL IS AWGG\n');
    case 2 % MWGG
        fprintf('\nTHE BEST FITTING MODEL IS MWGG\n');
    case 3 % POISSON
        fprintf('\nTHE BEST FITTING MODEL IS POISSON\n');        
end


% show in a Matlab figure all the estimated (mu,sigma) points 
% and the three fitting models (for the three noise models)
picturewidth = 20;
hw_ratio = 0.7; %Set for the height

fig_1 = figure(1);
plot(im_means(:),im_stdvs(:),'bo'); % point cloud
xlabel('$\mu$');ylabel('$\sigma$');
axis([0,255,0,16]);
legend('estimated ($\mu$,$\sigma$) points');

set(findall(fig_1,'-property','FontSize'),'FontSize',21);
set(findall(fig_1,'-property','Interpreter'),'interpreter','latex');
set(findall(fig_1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(fig_1,'units','centimeters','Position',[3 3 picturewidth hw_ratio * picturewidth]);
pos = get(fig_1,'Position');
set( fig_1,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] );
print(fig_1,'pdf_figure','-dpdf','-vector','-fillpage');  


fig_2 = figure(2);
plot(im_means(:),im_stdvs(:),'bo'); % point cloud
hold all;
% AWGG best fitting model (previously estimated)
plot([0,255],[AWGG_q,AWGG_q],'lineWidth',1.5,'Color','r');

% MWGG best fitting model (previously estimated)
plot([0,255],[0,255*MWGG_m],'lineWidth',1.5,'Color','k');%line that goes through (0,0) and (255,MWGG_m*255)

% POIS model (it did not need to be estimated)
plot(0:255,sqrt(0:255),'lineWidth',1.5,'Color','m');

xlabel('$\mu$');ylabel('$\sigma$');
axis([0,255,0,16]);
legend('estimated ($\mu$,$\sigma$) points','AWGG best fitting model','MWGG best fitting model','POIS model','Location','southeast');
title('Inference of the Noise')

set(findall(fig_2,'-property','FontSize'),'FontSize',21);
set(findall(fig_2,'-property','Interpreter'),'interpreter','latex');
set(findall(fig_2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(fig_2,'units','centimeters','Position',[3 3 picturewidth hw_ratio * picturewidth]);
pos = get(fig_2,'Position');
set( fig_2,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] );
print(fig_2,'pdf_figure','-dpdf','-vector','-fillpage');  


fig_100 = figure(100);
plot([0,255],[AWGG_q,AWGG_q],'LineWidth',1.5);
xlabel('$\mu$','Interpreter','Latex');ylabel('$\sigma$','Interpreter','Latex');
axis([0,255,0,16]);
title('AWGG Noise','Interpreter','latex');

set(findall(fig_100,'-property','FontSize'),'FontSize',21);
set(findall(fig_100,'-property','Box'),'Box','off');
set(findall(fig_100,'-property','Interpreter'),'interpreter','latex');
set(findall(fig_100,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(fig_1,'units','centimeters','Position',[3 3 picturewidth hw_ratio * picturewidth]);
pos = get(fig_100,'Position');
set( fig_100,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] );
print(fig_100,'pdf_figure','-dpdf','-vector','-fillpage');  

fig_101 = figure(101);
plot([0,255],[0,255*MWGG_m],'LineWidth',1.5);
xlabel('$\mu$','Interpreter','Latex');ylabel('$\sigma$','Interpreter','Latex');
axis([0,255,0,16]);
title('MWGG Noise','Interpreter','latex');

set(findall(fig_101,'-property','FontSize'),'FontSize',21);
set(findall(fig_101,'-property','Box'),'Box','off');
set(findall(fig_101,'-property','Interpreter'),'interpreter','latex');
set(findall(fig_101,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(fig_101,'units','centimeters','Position',[3 3 picturewidth hw_ratio * picturewidth]);
pos = get(fig_101,'Position');
set( fig_101,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] );
print(fig_101,'pdf_figure','-dpdf','-vector','-fillpage');  

fig_102 = figure(102);
plot(0:255,sqrt(0:255),'LineWidth',1.5);
xlabel('$\mu$','Interpreter','Latex');ylabel('$\sigma$','Interpreter','Latex');
axis([0,255,0,16]);
title('Poisson Noise','Interpreter','Latex');

set(findall(fig_102,'-property','FontSize'),'FontSize',21);
set(findall(fig_102,'-property','Box'),'Box','off');
set(findall(fig_102,'-property','Interpreter'),'interpreter','latex');
set(findall(fig_102,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(fig_102,'units','centimeters','Position',[3 3 picturewidth hw_ratio * picturewidth]);
pos = get(fig_102,'Position');
set( fig_102,'PaperPositionMode','Auto','PaperUnits','centimeters','papersize',[pos(3),pos(4)] );
print(fig_102,'pdf_figure','-dpdf','-vector','-fillpage');  


% if the "best fitting" model is AWGG or MWGG, we now estimate
% the shape parameter beta of the GG noise distribution (see slides). 
% To do that, we compute normalized histograms of estimated 
% noise realizations and then visually compare them to some
% theoretical GG distributions with different values of beta
if ( RMSRS_sorted_inds(1) <= 2 )
    %
    noise_H_bins_n = 49;

    im_stdv_max           = max(im_stdvs(:));
    noise_H_edges_min     = -3 * im_stdv_max;
    noise_H_edges_max     = +3 * im_stdv_max;
    noise_H_edges_n       = noise_H_bins_n + 1;
    noise_H_edges         = linspace(noise_H_edges_min,noise_H_edges_max,noise_H_edges_n);
    noise_H_bin_centers   = 0.5 * ( noise_H_edges( 1:(noise_H_edges_n - 1 ) ) + ...
                                    noise_H_edges( 2:noise_H_edges_n ) );
    noise_H_bins_size     = noise_H_edges(2) - noise_H_edges(1); 

    noise_pdf_xs = linspace(noise_H_edges_min,noise_H_edges_max,2000);
    betas = [1, 1.5, 2];
    noise_sigma    = mean(im_stdvs(:)) ;
    noise_pdf_GG_1 = generalized_gauss_pdf_1D(0,noise_sigma,betas(1),noise_pdf_xs); % mu, sigma, beta e xs sono i punti.
    noise_pdf_GG_2 = generalized_gauss_pdf_1D(0,noise_sigma,betas(2),noise_pdf_xs);
    noise_pdf_GG_3 = generalized_gauss_pdf_1D(0,noise_sigma,betas(3),noise_pdf_xs);

    pixel_noise_values = zeros(h,w,ims_n);
    for pixel_i = 1 : h
        for pixel_j = 1 : w
            pixel_values       = ims_seq(pixel_i,pixel_j,:);
            pixel_noise_values(pixel_i,pixel_j,:) = pixel_values - im_means(pixel_i,pixel_j);
        end
    end
    pixel_noise_values = pixel_noise_values(:);

    fig_4 = figure(4);
    pixel_noise_H      = histcounts( pixel_noise_values , noise_H_edges );
    pixel_noise_H      = pixel_noise_H ( 1:noise_H_bins_n ) / sum( pixel_noise_H( 1:noise_H_bins_n ) );
    pixel_noise_H      = pixel_noise_H / noise_H_bins_size;

    bar(noise_H_bin_centers,pixel_noise_H,'b'); hold on;
    plot(noise_pdf_xs,noise_pdf_GG_1,'r','linewidth',2);
    plot(noise_pdf_xs,noise_pdf_GG_2,'k','linewidth',2);
    plot(noise_pdf_xs,noise_pdf_GG_3,'g','linewidth',2);
    axis([noise_H_edges_min,noise_H_edges_max,0,1.05*max([max(pixel_noise_H),max(noise_pdf_GG_1)])]);
    legend('hist','GG pdf, $\beta=1$','GG pdf, $\beta=1.5$','GG pdf, $\beta=2$','interpreter','latex'); 
    xlabel('x','Interpreter','latex'); ylabel('p(x)','interpreter','Latex');
    title('Estimation of the $\beta$ parameter','interpreter','latex');

    set(findall(fig_4,'-property','FontSize'),'FontSize',21);
end 
    








