function C = generateConstraintMatrix(M,N)
% function C = geraMatRestr(M,NBF)
% Generates the broadside constraints matrix described in Frost - 1972 - An algorithm for
% Linearly Constrained Adaptive Array Processing 
% eqs (7) and (8)
% Parameters
%   M - number of microphones
%   NBF - length of the filters

if (M~=1)
    C = zeros (M*N,N);
% gera cada coluna de C a partir da equacao 8 de Frost - 1972
    for i=1:N
        C((i-1)*M+1:i*M, i) = ones (M,1);
    end,
else
    C = eye(N);
end,
