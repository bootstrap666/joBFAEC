function [HLOS, H] = generatechannels (M, Nh, F, T0, TR60, fs, PLOS, PREF)
% function [HLOS, H] = generatechannels (M, Nh, F, T0, TR60, fs, PLOS, PREF)
% Generates an echo impulse response matrix.
% Parameters
%   M           -   Number of microphones 
%   Nh          -   Resulting echo channels length
%   F           -   Number of delay samples between adjacent microphones for
%                   endfire signals
%   T0          -   Delay from the loudspeaker to the center of the microphone
%                   array
%   TR60        -   Reverberation time
%   fs          -   sampling rate
%   PLOS        -   Line of Sight (LOS) component total power
%   PREF        -   Reflections component total power
%   HLOS        -   LOS component echo channel matrix
%   H           -   Reflextions component echo channel matrix
[atrasos, deltaangulos] = anglepartition(F);

[HLOS, defasagem] = generateLOS(M,Nh,atrasos, deltaangulos, T0, fs, PLOS);
[H, hproto] = generatereflex (M, F, T0, TR60, defasagem, fs, ...
    PREF, Nh);

