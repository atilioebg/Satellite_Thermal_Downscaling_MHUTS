function [T] = TES(RAD, LAMB, EMI, LD, LU, Tr)

% Calculates the temperature from reference channel method
% IN:
% RAD  = radiance;
% LAMB = wavelenght in micrometer (without cientific notation(x10^y));
% EMI  = emissivity;
% LD   = downwelling radiation;
% 
% OUT:
% T    =  temperature

%**************************************************************************
%**************************************************************************
% c = 2.9979*10^8; % speed of light m/s
% h = 6.6261*10^-34; % Planck's constant J*s
% k = 1.3807*10^-23; % Boltzmann's constant J/K
%**************************************************************************
%**************************************************************************

C1 = 1.191042e8;    % <=== first radiation constant (W/um^4 m^2 sr) 

C2 = 1.4387752e4;   % <=== second radiation constant (mu K)

    for i = 1 : size(RAD,1)
        for j = 1 : size(RAD,2)
            
            A = (RAD(i,j) - (1 - EMI)*LD)/EMI;
            
            T(i,1) = C2/(LAMB*log((C1/((LAMB^5)*A))+1)); 

    end;
end