% Variables initialization
M               = 2;    % Number of microphones
F               = 4;    % Ratio between temporal and spatial oversampling factors
T0              = 1E-3; % approximate delay from loudspeaker to microphone from LOS path
TR60            = 1E-2; % room reverberation time
fs              = 8E3;  % sampling rate
PLOS            = 0.1;  % LEM Power gain in LOS
PREF            = 0.5;  % Total reflections LEM power gain
Nh              = ceil((T0 + TR60)*fs);  % LEM plants total length
Nbf             = 16; % Beamformer filter length
Naec            = Nh + Nbf -1; % AEC filter length (<Nh + Nbf -1)
SNR             = 20; % Far-end to noise ratio in dB
pctmucrit       = 0.05; % percentage of mucrit to use
realizations    = 10; % Number of realizations in Monte-Carlo simulation
iterations      = 1E5; % Number of iterations (input samples)
fCal            = [1; zeros(Nbf-1,1)];  % beamformer impulse response at the local 
                                        % speaker direction (broadside).
%%
disp('Loudspeaker-Enclosure-Microphone plant generation');
[HLOS, H] = generatechannels (M, Nh, F, T0, TR60, fs, PLOS, PREF);

%%
disp('step-size estimation');
% synthetic input signal parameters
% u[n] = alpha[n] + rhoU * u[n-1]
sigmaU = 1;
rhoU = -0.9;
% critical step-size computation (mucrit).
Rbb = generateRbb(H+HLOS, Naec, Nbf, sigmaU, SNR, rhoU);
C = generateConstraintMatrix(M,Nbf);
Ce = [zeros(Naec,Nbf); C];
Pe = eye (Naec+M*Nbf) - Ce*((Ce'*Ce)\Ce');
mucrit = 2/(3*trace(Pe*Rbb*Pe));
mu = 0.05*mucrit;
%%
disp('Initializing adaptive filter');
w0 = C*((C'*C)\fCal);
h0 = zeros(Naec,1);
%%
disp('Predicting the behavior with statistical model');
Dstatmodel = secondordermodeldiag(iterations, H+HLOS, mu, w0, h0, SNR, rhoU, fCal);
Jinf = calcJinf (H+HLOS, Naec, sigmaU, SNR, rhoU, fCal, mu);
%
fprintf('Monte Carlo Simulation (%d runs)\n',realizations);
[A, DMonteCarlo] = runjoBFAEC(realizations, iterations, H+HLOS, mu, mu, w0, h0, SNR, rhoU, ...
                    false, false, 'jobfaecdemo.mat');
%%                
disp('Comparing prediction and MC-Simulation');
plotsimresults(DMonteCarlo, Dstatmodel, Jinf)
