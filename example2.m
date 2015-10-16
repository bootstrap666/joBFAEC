function [M, NAEC, MU] = example2(H, mics, SNR, sigmaU, rhoU, fCal, Naecs, JdB1, niter1, pctmucrit)
% function [M, NAEC, MU] = example2(H, mics, SNR, sigmaU, rhoU, fCal, Naecs, JdB1, niter1, pctmucrit)
% Computes the MOP statistical model curves and designs all filters 
% with M(i) microphones that satisfies a MOP of Jdb1 in  niter1 iterations
% Parameters
%   H           -   Echo channel matrix 
%   mics        -   Candidate number of microphones vector
%   SNR         -   Echo to noise ratio in dB 
%   sigmaU      -   far-end signal variance
%   rhoU        -   far-end signal correlation 
%   fCal        -   impulse response in the desired direction 
%   Naecs       -   candidate AEC filter lengths 
%   JdB1        -   desired Mean Output Residual Echo Power at iteration nite1 
%   niter1      -   iteration at which an MOP of JdB1 is desired 
%   pctmucrit   -   percentage of the maximum step-size to be used
%   M           -   vector with the possible values of number of microphones 
%   NAEC        -   vector with the possible values of AEC filter length 
%   MU          -   vector with the possible values of step-size
Nh = size(H,1);

sigmaR = sigmaU*10^(-SNR/10);

Nbf = length(fCal);
J1 = 10^(JdB1/10);
M = [];
NAEC = [];
MU = [];
Ruhuh = sigmaU*toeplitz(rhoU.^(0:(Nh-1)),rhoU.^(0:(Nh-1)));


for i=1:length(mics)
    C = geraMatRestr(mics(i),Nbf);
    
    for j=1:length(Naecs)
        fprintf('Testing %d microphones with Naec=%d\n', mics(i),Naecs(j));
        if (mics(i) ~=1)
            Rbb = generateRbb(H(:,1:mics(i)), Naecs(j), Nbf, sigmaU, SNR, rhoU);

            Ce = [zeros(Naecs(j),Nbf); C];
        
            Omega = eye (Naecs(j)+mics(i)*Nbf) - Ce*((Ce'*Ce)\Ce');
            trac = trace(Omega*Rbb);
            % 5 percent of mucrit
            mu = pctmucrit * 2/(3*trac);
     
            aopt = Rbb\(Ce*((Ce'*(Rbb\Ce))\fCal));
        
            v = [zeros(Naecs(j),1); C*((C'*C)\fCal)];
            v = v - aopt;
            K = v*v';
        
            [Q, Lambda] = eigs(Omega*Rbb*Omega, Naecs(j) + (mics(i)-1)*Nbf);

            Jmin = aopt'*Rbb*aopt;

        else
                Rbb = sigmaU*toeplitz(sigmaU*(rhoU.^(0:(Naecs(j)-1))));
                Ruuh = sigmaU*toeplitz(rhoU.^(0:(Naecs(j)-1)),rhoU.^(0:(Nh-1)));
                
            trac = trace(Rbb);
            % 5 percent of mucrit
            mu = pctmucrit * 2/(3*trac);
            aopt = (Rbb\Ruuh) * H(:,1);
            K = (- aopt)*(- aopt)';            
            [Q, Lambda] = eig(Rbb);

            Jmin = sigmaR + H(:,1)'*Ruhuh*H(:,1) - H(:,1)'* Ruuh'*(Rbb\Ruuh) * H(:,1);
        end

        lambda = diag(Lambda);
        
        rhovv = diag(Q'*K*Q);
        
        B = diag((1-mu*lambda).^2 + (mu*lambda).^2) + mu^2*(lambda*lambda');
        
        rhovvinfty = mu^2*(eye (Naecs(j)+mics(i)*Nbf) - (B))\lambda;

        rhovv = (B^niter1) * (rhovv - rhovvinfty) + rhovvinfty;
        
        if ((rhovv'*lambda + Jmin) < J1)

            M = [M; mics(i)];
            NAEC = [NAEC; Naecs(j)];
            MU = [MU; mu];
        end,

    end,
end,

   
        