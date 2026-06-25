function [snr] = compute_snr(sig,ref)
%
% Compute Signal-to-Noise Ratio (SNR) defined by:
%
% snr(sig,ref) = 10 ln( num / den )   where:
%
% num = ||ref - ref_mean||_F^2
% den = ||ref - sig||_F^2
%
% and where || . ||_F denotes the Frobenious norm
%
% Usage:
%       snr = compute_snr(ref,sig)
%
% Input:
%       sig:  corrupted signal
%       ref:  reference signal
%  
% Output:
%       snr:  SNR value

num     = sum( ( ref(:) - mean(ref(:)) ).^2 );
den     = sum( ( ref(:) - sig(:)       ).^2 );
snr     = 10 * log10( num / den );

end


