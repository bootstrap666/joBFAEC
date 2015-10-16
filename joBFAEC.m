function aNew = joBFAEC(a, b, d, aq, Omega, MU, Naec, M, Nbf)
% function aNew = joBFAEC(a, b, d, aq, Omega, MU)
% Jointly optimized Beamformer-Acoustic echo Canceller weight update 
% equation
% anew - new weight vector
% a - current weight vector
% b - input vector 
%       b = [-u;  x ] where u is the far-end signal and x is the
%       microphone signals grouped according to Frost III (1976).
%       Function updateSampleVector can be used to 
% d - error signal
% aq - quiescent adjoint vector 
%       aq = [zeros (Naec,1); wq]
% Omega - Adjoint Projection matrix 
%       Omega = bdiag(eye(Naec), P), where P is the projection matrix in
%       the vector space orthogonal to the columns of constraint matrix C
% MU - step size matrix 
%       MU = diag ([muAEC*ones(1:Naec); muBF*ones(1:(M*Nbf))])
if (nargin < 7)
    aNew = Omega *(a -d*MU*b ) + aq;
else
    hnew = a(1:Naec) - d*MU(1,1)*b(1:Naec);
    P = Omega((Naec+1):end,(Naec+1):end);
    bnew = P *(a((Naec+1):end) -d*MU(end,end)*b((Naec+1):end) ) + aq((Naec+1):end);
    aNew = [hnew; bnew];
end