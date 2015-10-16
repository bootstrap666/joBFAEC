function fig = perfcurve1(H, mics, SNR, sigmaU, rhoU, points, fCal)
% function c1 = perfcurve1(H, mics, SNR, sigmaU, rhoU, points, fCal)
% Plots the steady-state MOP as a function of the step-size parameter.
% Assumes an AR(1) far-end signal and mu = mucrit
% Parameters
%   fig                 -   figure handle 
%   H                   -   LEM plant impulse response matrix (Nh x M)
%   mics                -   number of microphones vector to be tested 
%   SNR                 -   echo to noise ratio (in dB)
%   sigmaU              -   Far-end signal variance
%   rhoU                -   Far-end signal correlation 
%   points              -   number of step-sizes to test
%   fCal                -   impulse response in the look direction (assumed at the broadside)


sigmaR = sigmaU*10^(-SNR/10);
Nh = size(H,1);
Nbf = length(fCal);
fig = figure();
set(gca,'FontSize',14);
ColOrd = get(gca,'ColorOrder');

[m,n] = size(ColOrd);
xlabel('\mu','FontSize',14);
ylabel('E\{d^2[\infty]\}_{dB}','FontSize',14);
hold on;
M = [];
linstyle = {'-','--',':','-.'};

for i=1:length(mics)
    if (mics(i) ==1)
        if (rhoU == 0)
            Ruu = sigmaU*eye(Nh);
        else
            Ruu = sigmaU * toeplitz(rhoU.^(0:(Nh-1)));
        end,
        Jmin = sigmaR;
        trac = trace(Ruu);
        mucrit = 2/(3*trac);
        mu = mucrit*(1:points)/points;
        c1 = Jmin *(1+ (mu*trac)./(2-2*mu*trac));
    else
        Rbb = generateRbb(H(:,1:mics(i)), Nh+Nbf-1, Nbf, sigmaU, SNR, rhoU);
        aopt = coptsol (H(:,1:mics(i)), Nh + Nbf -1, sigmaU, SNR, rhoU, fCal);

        C = geraMatRestr(mics(i),Nbf);

        Ce = [zeros(Nh+Nbf-1,Nbf); C];
        
        Jmin = aopt'*Rbb*aopt;
        Omega = eye (length(aopt)) - Ce*((Ce'*Ce)\Ce');
        trac = trace(Omega*Rbb);
        mucrit = 2/(3*trac);
        mu = mucrit*(1:points)/points;
        c1 = Jmin *(1+ (mu*trac)./(2-2*mu*trac));
    
    end
    M = [M; sprintf('M = %d',mics(i))];
    Col = ColOrd(i,:);
    plot(mu, 10*log10(c1), char(linstyle(i)),'Color', Col, 'LineWidth', 2);
end
legend(M,'Location','Best');
grid on;
