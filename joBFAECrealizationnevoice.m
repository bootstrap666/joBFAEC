function [A, D] = joBFAECrealizationnevoice(iterations, systemTransferMatrix, ...
    aq, Omega, MU, decimationFactorC, Naec, M, Nbf, Nh, sigmaR, a0, ...
    showprogressbars, input)
% function [A, D] = joBFAECrealizationnevoice(iterations, systemTransferMatrix, ...
%     aq, Omega, MU, decimationFactorC, Naec, M, Nbf, Nh, sigmaR, a0, ...
%     showprogressbars, input)
% Monte Carlo simulation of a BF-AEC jointly optimized system with nonstationary 
% inputs due to the far-end input signal (single realization)
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
%   sigmaR                  -   local noise variance
%   a0                      -   extended weight vector initial coefficients
%   showprogressbars        -   show progress bars showing the progress of
%                               each iteration
%   voicevector             -   voice vector with at least iterations samples



A = zeros (Naec+M*Nbf, iterations/decimationFactorC + 1);
D = zeros (1, iterations);

c = zeros (Nh+(M+1)*Nbf -1,1);
c(Nh:-1:1) = input(length(input):-1:(length(input)-Nh+1));

a = a0;

A(:,1) = a0;

if (showprogressbars ==1)
    wb = waitbar(0, 'Delicious. Delicious. This way you''re gonna kill me');
end,

for i=1:iterations
% Far-End signal generation as an AR(1) process with correlation
% coefficient rhoU and variance sigmaU
    newFarEnd = input(i);

% Noise signal generation as a WGN with variance sigmaR
    newNoise = sqrt(sigmaR)*randn(M,1);
    
    [c, a, d, b] = joBFAECiterationt(newFarEnd, newNoise, ...
        systemTransferMatrix, c, a, aq, Omega, MU, Naec, Nbf, M,Nh);
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