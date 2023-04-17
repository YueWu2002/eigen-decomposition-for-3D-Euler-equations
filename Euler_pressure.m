function [p] = Euler_pressure(uu)
% compute the pressure
% 
% input: 
%   uu: conservative variables, size = [5, 1]
% 
% output: 
%   p:  presssure, positive scalar

gamma = 1.4;
gammahat = 0.4;

p = gammahat * (uu(5) - 0.5*sum(uu(2:4).^2) / uu(1));

if p <= 0.0
    error('Negative pressure detected!');
end

end