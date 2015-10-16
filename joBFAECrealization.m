function [A, D] = joBFAECrealization(iterations, systemTransferMatrix, ...
    aq, Omega, MU, decimationFactorC, Naec, M, Nbf, Nh, ...
    rhoU, sigmaU, sigmaR, a0, showprogressbars)
%function [A, D, zica] = joBFAECrealization(iterations, systemTransferMatrix, ...
%    aq, Omega, MU, decimationFactorC, Naec, M, Nbf, Nh, ...
%    rhoU, sigmaU, sigmaR, a0, showprogressbars)
% Monte Carlo simulation of a BF-AEC jointly optimized system (single realization)
% Parameters
%   A                       -   Weight vector regressor
%   D                       -   Instantaneous output power 
%   iterations              -   number of samples in each run of the algorithm
%   systemTransferMatrix    -   System Transfer Matrix as Described in
%                               Maruo, Resende and Bermudez 2013
%   aq                      -   extended quiescent vector
%   Omega                   -   extended Projection matrix 
%   MU                      -   step-size matrix 
%   decimationFactorC       -   Weight error matrix decimation factor (needed for long 
%                               LEM impulse responses to save memory)
%   Naec                    -   AEC Filter length
%   M                       -   Number of microphones 
%   Nbf                     -   Beamformer Filter length
%   Nh                      -   Echo channel impulse response length
%   rhoU                    -   far-end correlation coefficient
%   sigmaU                  -   far-end signal variance
%   sigmaR                  -   local noise variance
%   a0                      -   extended weight vector initial coefficients
%   showprogressbars        -   show progress bars showing the progress of
%                               each iteration
A = zeros (Naec+M*Nbf, iterations/decimationFactorC + 1);
D = zeros (1, iterations);

aU = [1; rhoU];

newFarEnd = 0;

c = zeros (Nh+(M+1)*Nbf -1,1);
%a = zeros(Naec+M*Nbf,1);

a = a0;

A(:,1) = a0;

if (showprogressbars ==1)
    wb = waitbar(0, 'Delicious. Delicious. This way you''re gonna kill me');
end,

%d2velho = 1E100;

% Fill the input vector 
oldNoise = zeros(M,1);
rhoR = 0;

for i=1:(Nh+Nbf-1)
    oldFarEnd = newFarEnd;
    newFarEnd = randn(1)*sqrt(1-rhoU^2);
    newFarEnd = sqrt(sigmaU)*(aU'*[newFarEnd; oldFarEnd]);

    newNoise = sqrt(sigmaR)*(randn(M,1) + rhoR*oldNoise);
    oldNoise = newNoise;
    c = updateSampleVectorc(c, newFarEnd, newNoise, Nh, Nbf, M);
end,
    
for i=1:iterations
% Far-End signal generation as an AR(1) process with correlation
% coefficient rhoU and variance sigmaU
    oldFarEnd = newFarEnd;
    newFarEnd = randn(1)*sqrt(1-rhoU^2);
    newFarEnd = sqrt(sigmaU)*(aU'*[newFarEnd; oldFarEnd]);
    %newFarEnd = beta* (1+sin(omega0*i))*newFarEnd; I dont remember why I
    %wrote this

% Noise signal generation as a WGN with variance sigmaR
    newNoise = sqrt(sigmaR)*(randn(M,1) + rhoR*oldNoise);
    oldNoise = newNoise;
%    newNoise = sqrt(sigmaR)*randn(M,1);
    [c, a, d, b] = joBFAECiteration(newFarEnd, newNoise, ...
        systemTransferMatrix, c, a, aq, Omega, MU, Nbf, M, Nh);
    if (~mod(i,decimationFactorC))
        A(:,i/decimationFactorC+1) = a;
        try 
            if (showprogressbars ==1)
                wb = waitbar(i/iterations,wb);
            end,
        end,
    end,
    
    D(i) = d^2;
end,
try
    if (showprogressbars==1)
        close (wb);
    end,
end,