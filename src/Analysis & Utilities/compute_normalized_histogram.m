function [H,H_bins_centers] = compute_normalized_histogram(xs,x_min,x_max,n_bins)

n_edges      = n_bins + 1;
edges        = linspace(x_min,x_max,n_edges);
bins_size    = edges(2) - edges(1); 

H = histcounts( xs , edges ); % absolute frequencies
H = H(1:n_bins) / sum( H(1:n_bins) ); % relative frequencies
H = H / bins_size; % normalized frequencies

H_bins_centers = 0.5 * ( edges(1:(n_edges-1)) + edges(2:n_edges) );

end

