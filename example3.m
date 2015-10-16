function example3(H, mics, Naecmax, SNR, sigmaU, rhoU, JinfdesdB, pctmucrit, fCal)
%function example3(H, mics, Naecmax, SNR, sigmaU, rhoU, JinfdesdB, pctmucrit, fCal)
% Computes the steady-state residual MOP as a function of the AEC length
% and prints the feasible configurations in the command window
% Parameters
%   H           -   Echo channel matrix 
%   mics        -   candidate microphone number vector
%   Naecmax     -   Maximum AEC length to test
%   SNR         -   Echo to noise ratio in dB 
%   sigmaU      -   far-end signal variance
%   rhoU        -   far-end signal correlation
%   JinfdesdB   -   Desired steady-state MOP in dB
%   pctmucrit   -   percentage of the maximum step-size to be used
%   fCal        -   Impulse response in the desired direction 

Nh = size(H,1);

Jinfdes = 10^(JinfdesdB/10);

sigmaR = sigmaU*10^(-SNR/10);

for i=1:length(mics)
    if (mics(i) == 1)
        Ruu = sigmaU * toeplitz(rhoU.^(0:(Naecmax-1)));
        Ruuh = sigmaU * toeplitz(rhoU.^(0:(Naecmax-1)), rhoU.^(0:(Nh-1)));
        Ruhuh = sigmaU * toeplitz(rhoU.^(0:(Nh-1)), rhoU.^(0:(Nh-1)));
        Jmin = sigmaR + H(:,1)'*Ruhuh*H(:,1) - H(:,1)'* Ruuh'*(Ruu\Ruuh) * H(:,1);
        if (Jmin > Jinfdes)
            continue;
        end,
        mu = 2*(Jinfdes - Jmin)/(Jinfdes*trace(Ruu));
        
        if (mu > pctmucrit*2/(3*trace(Ruu)))
            mu = pctmucrit*2/(3*trace(Ruu));
        end,
        
        Jinf =  Jmin * mu *trace(Ruu) / (2 - mu *trace(Ruu));
        fprintf('M=%d, mu=%e, Jinf(dB) = %f\n', mics(i), mu, 10*log10(Jinf + Jmin));
    else
        C = geraMatRestr(mics(i),length(fCal));
        Ce = [zeros(Naecmax,length(fCal)); C];
        Omega = eye (Naecmax + mics(i)*length(fCal)) - Ce*((Ce'*Ce)\Ce');
        Rbb = generateRbb(H(:,1:mics(i)), Naecmax, length(fCal), sigmaU, SNR, rhoU);
        aopt = coptsol (H(:,1:mics(i)), Naecmax, sigmaU, SNR, rhoU, fCal);
        Jmin = aopt'*Rbb*aopt;
        if (Jmin > Jinfdes)
            continue;
        end,
        trac = trace(Omega*Rbb*Omega);
        mu = 2*(Jinfdes - Jmin)/(Jinfdes*trac);
        
        if (mu > pctmucrit*2/(3*trac))
            mu = pctmucrit*2/(3*trac);
        end,
        
        Jinf =  Jmin * mu *trac / (2 - mu *trac);
        fprintf('M=%d, mu=%e, Jinf(dB) = %f\n', mics(i), mu, 10*log10(Jinf + Jmin));
        
    end,

end,

