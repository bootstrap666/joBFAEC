function [A, D] = runjoBFAEC(realizations, iterations, H, muAEC, muBF, w0, h0, SNR, rhoU, showprogressbars, continuefrom, filename)
% function [A, D] = runjoBFAEC(realizations, iterations, H, muAEC, muBF, 
% w0, h0, SNR, rhoU, showprogressbars, continuefrom, filename)
% Monte Carlo simulation of a BF-AEC jointly optimized system
% Parameters
%   A                   -   Weight vector regressor
%   D                   -   Mean output power 
%   realizations        -   number of runs of the algorithm
%   iterations          -   number of samples in each run of the algorithm
%   H                   -   LEM plant impulse response matrix (Nh x M)
%   muAEC               -   AEC step-size
%   muBF                -   beamformer step-size
%   w0                  -   beamformer initial coefficients
%   h0                  -   AEC initial coefficients
%   SNR                 -   echo to noise ratio
%   rhoU                -   far-end correlation coefficient
%   showprogressbars    -   show progress bars showing the progress of each 
%                           algorithm run
%   continuefrom        -   load the simulation state from a file 
%   filename            -   simulation state file

sigmaU = 1;
% convert the noise power from dB scale to linear scale
sigmaR = sigmaU*10^(-SNR/10);

% Number of microphones/length of the LEM plants
[Nh, M] = size(H);
% AEC length
Naec = length(h0);
% beamformer filter length
Nbf = length(w0)/M;

% decimation factor to save memory
decimationFactorC = iterations / 100;
%decimationFactorD = iterations / 100;

% weight vector initialization
a0 = zeros(Naec + M*Nbf,1);
a0(1:Naec) = h0;
a0((Naec+1):(Naec+M*Nbf)) = w0;


% Quiescent adjoint vector 
C = generateConstraintMatrix(M,Nbf);
% Desired response in the DOA (0 rad)
fCal = [1; zeros(Nbf-1,1)];

%Example using Buckley's eigenvector constraints
%angle = pi/4;
%desiredimpres = sinc((-(Nbf-1)/2):((Nbf-1)/2));
%dmic = 5E-3;
%c = 343;
%fs = 2^13;
%threshold = 0.99;

%[C, fCal] = generateRestrMatrixAng (angle, desiredimpres, dmic, c,...
%     fs, M, threshold);

% quiescent vector
wq = C*((C'*C)\fCal);

% initial weights
aq = [zeros(Naec,1);wq];

% System Transfer Matrix
systemTransferMatrix = transferMatrix(H, Naec, Nbf);

% Step size Matrix

MU = diag([muAEC*ones(Naec,1); muBF*ones(M*Nbf,1);]);

% Adjoint projection Matrix
% Omega = zeros(Naec+M*Nbf, Naec + M*Nbf);
% Omega (1:Naec,1:Naec) = eye(Naec);
% P = eye (M*Nbf)-C*((C'*C)\C');
% Omega ((Naec+1):(Naec+M*Nbf),(Naec+1):(Naec+M*Nbf)) = P;
Ce = [zeros(Naec,Nbf);C];
Omega = eye(length(aq)) - Ce*((Ce'*Ce)\Ce');

A = zeros (length(aq),iterations/decimationFactorC+1);
D = zeros (1,iterations);

if (~continuefrom)
    i=1;
else
    load (filename);
end,
while (i<=realizations)
    fprintf('realization %d\n',i);
    [Ar, Dr] = joBFAECrealization(iterations, ...
        systemTransferMatrix, aq, Omega, MU, decimationFactorC, Naec, ...
        M, Nbf, Nh, rhoU, sigmaU, sigmaR, a0, showprogressbars);
    A = A + Ar;
    D = D + Dr;
    i = i+1;
    save (filename);
end,
A = A/realizations;
D = D/realizations;