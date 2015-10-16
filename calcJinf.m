function Jinf = calcJinf (H, Naec, sigmaU, SNR, rhoU, fCal, mu)
% function Jinf = calcJinf (H, Naec, Nbf, sigmaU, SNR, rhoU, fCal, mu)
% Estimates the steady-state residual error of a jointly optimized BF-AEC
% echo cancellation system with the LMS-like adaptation.
% The model assumes the far-end signal is an AR(1) process and the desired
% DOA is orthogonal to an Uniform linear array of microphones.
% H         - LEM plants impulse response matrix
% Naec      - Acoustic echo canceller length 
% sigmaU    - far-end input power
% SNR       - far-end to noise ratio (dB)
% rhoU      - far-end autocorrelation coefficient
% fCal      - desired impulse response in the direction orthogonal to the
%               array
% mu        - algorithm step-size

Nbf =  length (fCal);
M = size(H,2);

C = generateConstraintMatrix(M,Nbf);
% Example using Buckley's Constraint matrix instead
% angle = pi/4;
% desiredimpres = sinc((-(Nbf-1)/2):((Nbf-1)/2));
% dmic = 5E-3;
% c = 343;
% fs = 2^13;
% threshold = 0.99;
% 
% [C, fCal] = generateRestrMatrixAng (angle, desiredimpres, dmic, c,...
%      fs, M, threshold);

%Extended constraint matrix
Ce = [zeros(Naec,Nbf); C];

Rbb = generateRbb(H, Naec, Nbf, sigmaU, SNR, rhoU);

aopt = Rbb\(Ce*((Ce'*(Rbb\Ce))\fCal));

Jmin = aopt'*Rbb*aopt;

Omega = eye (length(aopt)) - Ce*((Ce'*Ce)\Ce');

% Excess error calculated from Maruo,Resende,Bermudez, 2013,Statistical 
% Analysis of a Jointly Optimized Beamformer-Assisted Acoustic Echo Canceler, 
% IEEE Transactions in Signal Processing
Jex = Jmin*mu*trace(Omega*Rbb*Omega) / (2 - mu*trace(Omega*Rbb*Omega));

Jinf = Jex + Jmin;

