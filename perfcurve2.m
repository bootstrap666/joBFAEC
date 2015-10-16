function P2 = perfcurve2(H, mics, SNR, sigmaU, rhoU, Naecmin, fCal, pctmucrit)
% function P2 = perfcurve2(H, mics, SNR, sigmaU, rhoU, Naecmin, fCal, pctmucrit)
% Plots the steady-state MOP as a function of the AEC length.
% Assumes an AR(1) far-end signal and mu = pctmucrit*mucrit
% Parameters
%   P2          -   MOP matrix  
%   H           -   LEM plant impulse response matrix (Nh x M)
%   mics        -   number of microphones vector to be tested 
%   SNR         -   echo to noise ratio (in dB)
%   sigmaU      -   Far-end signal variance
%   rhoU        -   Far-end signal correlation 
%   Naecmin     -   minimum AEC length to test
%   fCal        -   impulse response in the look direction (assumed at the broadside)
%   pctmucrit   -   percentage of the critical step-size

sigmaR = sigmaU*10^(-SNR/10);
Nh = size(H,1);
Nbf = length(fCal);
figure();
set(gca,'FontSize',14);
ColOrd = get(gca,'ColorOrder');

[m,n] = size(ColOrd);
xlabel('N_{AEC}','FontSize',14);
ylabel('E\{d^2[\infty]\}_{dB}','FontSize',14);
hold on;
M = [];
h0 = H(:,1);
xlim([Naecmin, Nh])
if (rhoU == 0)
    Ruhuh = sigmaU*eye(Nh);
else
    Ruhuh = sigmaU * toeplitz(rhoU.^(0:(Nh-1)));
end

P2 = [];
NaecMax = (Nh+Nbf-1);
%NaecMax = 1000;

for k=1:length(mics)
    for i=Naecmin:NaecMax
        fprintf('M = %d, Naec = %d\n',mics(k),i);
        if (mics(k)==1)
            if (rhoU == 0)
                Ruu = sigmaU*eye(i);
                Ruuh = [Ruu, zeros(i, Nh - i)];
            else
                Ruu = sigmaU * toeplitz(rhoU.^(0:(i-1)));
                Ruuh = sigmaU * toeplitz(rhoU.^(0:(i-1)),rhoU.^(0:(Nh-1)));
            end,
            Jmin = sigmaR + h0'*Ruhuh*h0 - h0'* Ruuh'*(Ruu\Ruuh) * h0;
            trac = trace(Ruu);
            mu = pctmucrit*2/(3*trac);
            c1(i-Naecmin+1) = Jmin*(1+(mu*trac)/(2- mu*trac));
        else
            Rbb = generateRbb(H(:,1:mics(k)), i, Nbf, sigmaU, SNR, rhoU);
%            aopt = coptsol (H(:,1:mics(k)), i, sigmaU, SNR, rhoU, fCal);

            C = geraMatRestr(mics(k),Nbf);

            Ce = [zeros(i,Nbf); C];
            
            aopt = Rbb\(Ce*((Ce'*(Rbb\Ce))\fCal));
        
            Jmin = aopt'*Rbb*aopt;
            Omega = eye (length(aopt)) - Ce*((Ce'*Ce)\Ce');
            trac = trace(Omega*Rbb);
            mu = pctmucrit*2/(3*trac);
            c1(i-Naecmin+1) = Jmin *(1+ (mu*trac)./(2-mu*trac));
        end,
    end
    P2 = [P2; c1];
    M = [M; sprintf('M = %d',mics(k))];
    Col = ColOrd(k,:);
    % Plot the data
    plot(Naecmin:NaecMax, 10*log10(c1),'Color', Col, 'LineWidth', 2);
end
legend(M,'Location','Best');
grid on;

        
  