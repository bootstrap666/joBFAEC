function [M, Naec, mu, Jinf] = example2p2(H, Mposs, SNR, sigmaU, rhoU, fCal, NAECS,  MU)
% function [M, Naec, mu, Jinf] = example2p2(H, Mposs, SNR, sigmaU, rhoU, fCal, NAECS,  MU)
% Uses the results from example2.m to select the filters with the minimal
% steady-state MOP for each number of microphones in Mposs
% Parameters
%   H           -   Echo channel matrix 
%   Mposs       -   candidate microphone number vector
%   SNR         -   Echo to noise ratio in dB 
%   sigmaU      -   far-end signal variance
%   rhoU        -   far-end signal correlation
%   fCal        -   Impulse response in the desired direction 
%   NAECS       -   Candidate AEC lengths
%   MU          -   Candidate Step-sizes
%   M           -   Resulting number of microphones vector
%   Naec        -   Resulting AEC length vector corresponding to M
%   mu          -   Resulting Step-Sizes vector corresponding to M 
%   Jinf        -   Resulting steady-state residual echo MOP vector corresponding
%                   to M
Jinf = 1E99;
Nbf = length(fCal);

for i=1:length(Mposs)
    if (Mposs(i) ~=1)
        Jinftemp = calcJinf (H(:,1:Mposs(i)), NAECS(i), Nbf, sigmaU, SNR, rhoU, fCal, MU(i));
    else
        Ruhuh = sigmaU*toeplitz(rhoU.^(0:(Nh-1)),rhoU.^(0:(Nh-1)));
        Ruu = sigmaU*toeplitz(sigmaU*(rhoU.^(0:(NAECS(i)-1))));
        Ruuh = sigmaU*toeplitz(rhoU.^(0:(NAECS(i)-1)),rhoU.^(0:(Nh-1)));
        Jmin = sigmaR + H(:,1)'*Ruhuh*H(:,1) - H(:,1)'* Ruuh'*(Ruu\Ruuh) * H(:,1);
        trac = trace(Ruu);
        Jinftemp = Jmin* (1 + MU(i)*trac/(2 - MU(i)*trac));
    end,
    if (Jinftemp < Jinf)
        M = Mposs(i);
        Naec = NAECS(i);
        mu = MU(i);
        Jinf = Jinftemp;
    end,
end