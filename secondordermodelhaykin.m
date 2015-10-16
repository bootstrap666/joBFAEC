function Dmoddiaglms = secondordermodelhaykin(iterations, H, mu, h0, SNR, rhoU)
% function Dmoddiag = secondordermodelhaykin(iterations, H, mu, h0, SNR, rhoU)
% Second moment statistics of the LMS algorithm in the system
% identification configuration with unknown plant impulse response H.
% Assumes an AR(1) stationary process as the far-end input.
%   Dmoddiag    -   Excess echo + noise + local speech output E{d^2[n]}
%                   statistics decimated by a factor of 100
%   iterations  -   simulation number of iterations 
%   H           -   LEM impulse response vector 
%   mu          -   adaptation step-size
%   h0          -   acoustic echo canceller initial weights
%   SNR         -   local speech + noise power to far-end power ratio 
%                       E{(local speech + noise)^2} / E{u^2[n]}
%   rhoU        -   AR(1) parameter -a_1. Controls the correlation of the
%                       far-end signal

sigmaU = 1;

%aopt = H(:,1);
Nh = size(H,1);
Naec = length(h0);

sigmaR = sigmaU*10^(-SNR/10);

%rho = sigmaU*rhoU.^[0:(Naec-1)];

decimationFactorC = iterations / 100;

if (iterations>100)
    decimationFactorD = 100;
else
    decimationFactorD = 1;
end,


%if (iterations>1E6)
%    decimationFactorD = iterations / 1E6;
%else
%    decimationFactorD = 1;
%end,
Ruu = sigmaU * toeplitz(rhoU.^(0:(Naec-1)));
Ruuh = sigmaU * toeplitz(rhoU.^(0:(Naec-1)),rhoU.^(0:(Nh-1)));
Ruhuh = sigmaU * toeplitz(rhoU.^(0:(Nh-1)),rhoU.^(0:(Nh-1)));
%Ruhuh = sigmaU*rhoU * toeplitz(rhoU.^(0:(Nh-1)));
aopt = (Ruu\Ruuh) * H(:,1);

Jmin = sigmaR + H(:,1)'*Ruhuh*H(:,1) - H(:,1)'* Ruuh'*(Ruu\Ruuh) * H(:,1);
%Jmin = sigmaR;

%v = h0; 
v = h0 - aopt;
K = v*v';

Dmoddiaglms = zeros( iterations/decimationFactorD ,1);
%MU = zeros (length(v));
%MU(1:Naec, 1:Naec) = muAEC * eye(Naec);
%MU((Naec+1):(Naec+M*Nbf), (Naec+1):(Naec+M*Nbf)) = muBF * eye(M*Nbf);

[Q, Lambda] = eig(Ruu);

lambda = diag(Lambda);

B = diag((1-mu*lambda).^2 + (mu*lambda).^2) + mu^2*(lambda*lambda');

rhovv = diag(Q'*K*Q);

try
    wb = waitbar(0,'Bring the caçamba. Bring the caçamba. Throw it all in there');
end,
Dmoddiaglms(1) = trace(K*Ruu) + Jmin;
% A = MU*Omega*Rbb;
% B = Omega*Rbb;
% C = Omega*MU;
% F = A*C;

%%
% A = MU*Rbb;
% B = Rbb;
% C = MU;
% F = A*C;
%%
%A = MU*Rbb;
%B = A*MU;

for i=2:iterations-1
%     K = K - A*K*Omega - B*K*C ...
%         + A*K*A' ...
%         + (trace(K*Rbb) + Jmin)*F;

    rhovv = B*rhovv + mu^2*Jmin*lambda;
    
    if (~mod(i,decimationFactorD))
        Dmoddiaglms(i/decimationFactorD+1) = rhovv'*lambda+Jmin;
    end,

    if (~mod(i,decimationFactorC))
        sprintf ('iteration %d\n',i);
        try 
            wb = waitbar(i/iterations,wb);
        end
    end,
end,
try 
    close (wb);
end,
