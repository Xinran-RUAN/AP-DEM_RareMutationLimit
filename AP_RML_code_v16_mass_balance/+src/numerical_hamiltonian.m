function Hhat = numerical_hamiltonian(u, dtheta, method, alpha)
%NUMERICAL_HAMILTONIAN Approximate -|u_theta|^2.
if nargin < 3 || isempty(method)
    method = 'godunov';
end
if nargin < 4 || isempty(alpha)
    alpha = 10;
end
u = u(:).';
switch lower(method)
    case 'godunov'
        [dm, dp] = src.one_sided_derivatives(u, dtheta, 'first-order');
        Hhat = -(min(dm,0)).^2 - (max(dp,0)).^2;
    case {'lf','lax-friedrichs','laxfriedrichs'}
        [dm, dp] = src.one_sided_derivatives(u, dtheta, 'first-order');
        Hhat = -((dm + dp)/2).^2 + 0.5*alpha*(dm - dp);
    case 'weno5'
        [dm, dp] = src.one_sided_derivatives(u, dtheta, 'weno5');
        Hhat = src.roe_flux_minus_p2(dp, dm);
    otherwise
        error('Unknown numerical Hamiltonian "%s".', method);
end
Hhat = Hhat(:).';
end
