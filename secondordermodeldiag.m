function Dmoddiag = secondordermodeldiag(iterations, H, mu, w0, h0, SNR, rhoU, fCal)
% function Dmoddiag = secondordermodeldiag(iterations, H, mu, w0, h0, ...
%                   SNR, rhoU, fCal)
% Second moment statistics of the jointly-optimized BF-AEC scheme based on
% the Frost Algorithm. Desired direction of arrival is assumed normal to
% the uniform linear array of microphones
%   Dmoddiag    - Excess echo + noise + local speech output E{d^2[n]}
%                   statistics decimated by a factor of 100
%   iterations  - simulation number of iterations 
%   H           - LEM impulse response matrix with each column containing 
%                   the LEM impulse response to one microphone
%   mu          - adaptation step-size
%   w0          - beamformer initial weights
%   h0          - acoustic echo canceller initial weights
%   SNR         - local speech + noise power to far-end power ratio 
%                   E{(local speech + noise)^2} / E{u^2[n]}
%   rhoU        - AR(1) parameter -a_1. Controls the correlation of the
%                   far-end signal
%   fCal        - desired impulse response in the desired direction of
%                   arrival

showwaitbar=0;
% variables initialization
sigmaU = 1;

M = size(H,2);
Naec = length(h0);
Nbf = length(w0)/M;

decimationFactorC = iterations / 100;
if (iterations>10000)
%if (iterations>100)
    decimationFactorD = 100;
else
    decimationFactorD = 1;
end,

C = generateConstraintMatrix(M,Nbf);

% For directions different than normal to the ULA the Buckley method to
% obtain eigenvalue constraints with minimal least squares error can be 
% used
%angle = pi/4;
% % Assuming a flat frequency response with Nbf/2 samples delay
%desiredimpres = sinc((-(Nbf-1)/2):((Nbf-1)/2));
%dmic = 5E-3;
%c = 343;
%fs = 2^13;
%threshold = 0.99;

%[C, fCal] = generateRestrMatrixAng (angle, desiredimpres, dmic, c,...
%     fs, M, threshold);

% The extended constraint matrix
Ce = [zeros(Naec,Nbf); C];

% Autocorrelation matrix Rss
Rbb = generateRbb(H, Naec, Nbf, sigmaU, SNR, rhoU);

% Optimal solution as derived on Frost III (1972)
aopt = Rbb\(Ce*((Ce'*(Rbb\Ce))\fCal));

% Optimal residual echo power + noise + local speech
Jmin = aopt'*Rbb*aopt;

% model parameters
v = [h0; w0];
v = v - aopt;
K = v*v';

Dmoddiag = zeros( iterations/decimationFactorD ,1);

Omega = eye (length(v)) - Ce*((Ce'*Ce)\Ce');

OmegaRbbOmega = Omega*Rbb*Omega;
OmegaRbbOmega = (OmegaRbbOmega + OmegaRbbOmega')/2;

%Jinf = Jmin*mu*trace(OmegaRbbOmega) (2 - trace(OmegaRbbOmega));

[Q, Lambda] = eigs(OmegaRbbOmega, Naec + (M-1)*Nbf);

lambda = diag(Lambda);

B = diag((1-mu*lambda).^2 + (mu*lambda).^2) + mu^2*(lambda*lambda');

rhovv = diag(Q'*K*Q);
if (showwaitbar)
    try
%         
        wb = waitbar(0,'Bring the caçamba. Bring the caçamba. Throw it all in there');
    end,
end,

Dmoddiag(1) = trace(K*OmegaRbbOmega) + Jmin;

%while ((10*log10(rhovv'*lambda+Jmin) - Jinf)>0.2)
for i=2:iterations-1

    rhovv = B*rhovv + mu^2*Jmin*lambda;
    
    if (~mod(i,decimationFactorD))
        Dmoddiag(i/decimationFactorD+1) = rhovv'*lambda+Jmin;
    end,

    if (~mod(i,decimationFactorC))
        sprintf ('iteration %d\n',i);
        if (showwaitbar)
            try 
                wb = waitbar(i/iterations,wb);
            end
        end,
    end,
end,
if (showwaitbar)
    try 
        close (wb);
    end,
end,
