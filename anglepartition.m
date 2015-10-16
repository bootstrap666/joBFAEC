function [delays, deltaangle] = anglepartition(F)
%function [delays, deltaangle] = anglepartition(F)
% Computes the echo power density of an uniform echo field as a function of
% the number of samples delay between adjacent microphones
% Parameters
%   F           -   number of delay samples between adjacent microphones
%                   for endfire signals
%   delays      -   vector of delay samples
%   deltaangle  -   echo power density

delays = [-1:-1:-F -(F-1):0 1:F (F-1):-1:0];

theta0 = asin((-F:(F-1))/F+1/(2*F));

angles = [-pi+theta0(F+1:2*F) theta0 pi-theta0(2*F:-1:F+1)];

deltaangle = [angles(2:4*F)-angles(1:4*F-1) 2*(pi-angles(4*F))];
