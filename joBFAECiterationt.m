function [c, a, d, b] = joBFAECiterationt(newFarEnd, newNoise, ...
    systemTransferMatrix, c, a, aq, Omega, MU, Naec, Nbf, M,Nh)
% function [snew, wnew, d, b] = joBFAECiterationt(newFarEnd, newNoise, ...
%     systemTransferMatrix, s, w, aq, Omega, MU, Naec, Nbf, M, Nh)
% Computes an iteration of the Jointly-Optimized Beamformer-Assisted Acoustic Echo Canceller 
% Parameters (faster version for diagonal block step-size matrix)
%   newFarEnd               -   New far-end signal sample 
%   newNoise                -   New beamformer noise vector (M x 1)
%   systemTransferMatrix    -   System Transfer Matrix as Described in
%                               Maruo, Resende and Bermudez 2013
%   ur                      -   Old input vector
%   w                       -   Old weight coefficients vector
%   aq                      -   extended quiescent vector
%   Omega                   -   extended Projection matrix 
%   MU                      -   step-size matrix 
%   Naec                    -   AEC Filter length
%   Nbf                     -   Beamformer Filter length
%   M                       -   Number of microphones 
%   Nh                      -   Echo channel impulse response length
%   urnew                   -   New input vector
%   wnew                    -   New weight coefficients vector
%   d                       -   Excess echo power
%   s                       -   New measured signal vector

c = updateSampleVectorc(c, newFarEnd, newNoise, Nh, Nbf, M);
b = systemTransferMatrix*c;

d = a'*b;

a = joBFAEC(a, b, d, aq, Omega, MU, Naec, M, Nbf);
