function Rbb = generateRbb(H, Naec, Nbf, sigmaU, SNR, rhoU)
% function Rbb = generateRbb(H, Naec, Nbf, sigmaU, SNR, rhoU)
% Generates the input signal analytical correlation matrix assuming an
% AR(1) far-end input
% Parameters
%   H           -   Echo channel matrix 
%   Naec        -   AEC impulse response length to test
%   Nbf         -   Beamformer filters length
%   sigmaU      -   far-end signal variance
%   sigmaU      -   far-end signal variance
%   SNR         -   Echo to noise ratio in dB 
%   rhoU        -   far-end signal correlation
%   Rbb         -   Input signal analytical correlation matrix

M = size(H,2);
Nh = size(H,1);
sigmaR = sigmaU*10^(-SNR/10);

systemTransferMatrix = transferMatrix(H, Naec, Nbf);

Rcc = zeros (Nh+Nbf + M*Nbf-1, Nh+Nbf + M*Nbf-1);

rho = sigmaU* (rhoU.^(0:(Nh+Nbf-2)));

Rcc(1:(Nh+Nbf-1),1:(Nh+Nbf-1)) = toeplitz(rho,rho);

Rcc((Nh+Nbf):(Nh+Nbf + M*Nbf-1),(Nh+Nbf):(Nh+Nbf + M*Nbf-1)) = ...
    sigmaR*eye(M*Nbf,M*Nbf);

Rbb = systemTransferMatrix*Rcc*systemTransferMatrix';

