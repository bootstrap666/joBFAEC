function cNew = updateSampleVectorc(c, newFarEnd, newNoise, Nh, Nbf, M)
% function bNew = updateSampleVector(c, newFarEnd, newNoise)
% Parameters
%   bNew        -    New sample vector
%   b           -   Current sample vector
%   newFarEnd   -   far-end signal current sample
%   newNoise    -   perturbation signal current sample
cNew = [newFarEnd; c(1:(Nh+Nbf-2)); newNoise; c((Nh+Nbf):(length(c)-M))];