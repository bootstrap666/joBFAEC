function [H, hproto] = generatereflex (M, F, T0, TR60, LOSdelay, ...
    fs, power, N)
%function [H, hproto] = generatereflex (M, F, T0, TR60, LOSdelay, ...
%    fs, power, N)
% Generate reflection component echoes 
% Assumes:
%   - Uniform echo power density across all arriving angles
%   - echo impulse responses are independent for different directions of
%   arrival
%   - Exponential model for echo impulse responses
% Parameters
%   M           -   Number of microphones
%   F           -   Number of delay samples between adjacent microphones for
%                   endfire signals
%   T0          -   Delay from the loudspeaker to the center of the microphone
%                   array
%   TR60        -   Reverberation time
%   LOSdelay    -   Delay from random LOS direction
%   fs          -   sampling rate
%   power       -   Total power gain from echo reflections
%   N           -   Echo impulse response length
if nargin < 6
    fs = 8192;
end,

if nargin < 7
    power = 1/2;
end,

if nargin < 8
    N = 2^nextpow2(ceil((TR60+T0)*fs));
end,

[atrasos, deltaangulos] = anglepartition(F);

deltaangulosnormalizado = deltaangulos/sum(deltaangulos);

n0 = ceil(T0*fs);

n0LOS = max (n0, n0-(M-1)*LOSdelay);

hproto = generateindepchannels (length(deltaangulos), n0, TR60, fs, ...
    power*deltaangulosnormalizado, N);


H = zeros(N,M);

for k=1:M
    n0k = n0LOS + (k-1)*LOSdelay;
    for i=1:length(deltaangulos)
        am_iniciohproto = n0k + (k-1)*atrasos(i);
        am_fimhproto = min(N, am_iniciohproto+(N-n0k));
        H(n0k:N,k) = H(n0k:N,k) + [zeros(am_iniciohproto-n0k,1); ...
            hproto(am_iniciohproto:am_fimhproto,i) ];% ...
    end,
end,
