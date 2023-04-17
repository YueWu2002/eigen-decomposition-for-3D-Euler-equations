function [Fn, d, L, R, J] = Euler_directional_flux(uu, nm)
% compute the characteristic decomposition of the flux Jacobian in a given direction
% 
% input: 
%   uu:     conservative variables, size = [5, 1]
%   nm:     non-zero vector that represents a direction, size = [3, 1] or [1, 3]
% 
% output: 
%   Fn:     directional flux, size = [5, 1]
%   J:      directional flux Jacobian, size = [5, 5]
%   d:      vector of eigen values, size = [5, 1]
%   L, R:   left and right eigen vectors, L = R^(-1), size = [5, 5]
% 
% CHECKED.

if nargout <= 0
    return;
end

gamma = 1.4;
gammahat = 0.4;

% compute a transform (represented by a 3*3 orthogonal matrix Q)
[Q, R] = qr(nm);
if R(1,1) < 0.0
    Q(:,1) = -Q(:,1);
end
nm = Q(:,1);

p = Euler_pressure(uu);
vel = uu(2:4) ./ uu(1);
H = (uu(5) + p) / uu(1);
c = sqrt(gamma * p / uu(1));

u2 = dot(vel, vel);
u2c = u2/c;
mach = vel ./ c;            % component-wise mach number

mn = dot(uu(2:4), nm(:));   % normal momentum
vn = mn / uu(1);            % normal velocity

% flux function in the given direction
Fn = [mn; p*nm(:) + vn*uu(2:4); H*mn];

if nargout >= 2
    % eigen values / wave speeds (left sonic wave, density contact discontinuity, tangential velocity contact discontinuity, right sonic wave)
    d = [vn - c; vn; vn; vn; vn + c];
end

if nargout >= 3
    % right and left eigen vectors (L*R = I)
    L = zeros(5,5);
    R = zeros(5,5);

    % left sonic wave
    R(:,1) = [1/c; mach(:) - nm(:); 0.5*u2c + c/gammahat - vn];
    L(1,:) = [0.25*gammahat*u2c + 0.5*vn; -0.5*gammahat*mach(:) - 0.5*nm(:); 0.5*gammahat/c];

    % right sonic wave
    R(:,5) = [1/c; mach(:) + nm(:); 0.5*u2c + c/gammahat + vn];
    L(5,:) = [0.25*gammahat*u2c - 0.5*vn; -0.5*gammahat*mach(:) + 0.5*nm(:); 0.5*gammahat/c];

    % density contact discontinuity
    R(:,2) = [1/c; mach(:); 0.5*u2c];
    L(2,:) = [c - 0.5*gammahat*u2c; gammahat*mach(:); -gammahat/c];

    % tangential velocity contact discontinuity
    vel_tilde = Q' * vel(:);
    R(:,3) = [0; Q(:,2); vel_tilde(2)];
    R(:,4) = [0; Q(:,3); vel_tilde(3)];
    L(3,:) = [-vel_tilde(2); Q(:,2); 0];
    L(4,:) = [-vel_tilde(3); Q(:,3); 0];
end


if nargout >= 5
    % flux Jacobian
    J = zeros(5, 5);
    J(:,1) = [0; 0.5*gammahat*u2*nm(:) - vn*vel(:); 0.5*gammahat*u2*vn - vn*H];
    J(:,2:4) = [nm(:)'; vn.*eye(3) + kron(vel(:), nm(:)') - gammahat*kron(nm(:), vel(:)'); (H*nm(:) - gammahat*vn*vel(:))'];
    J(:,5) = [0; gammahat*nm(:); gamma*vn];
end

end