function T = upwind_transport_matrix(u, op)
%UPWIND_TRANSPORT_MATRIX Matrix for D_theta F, F=-W u_theta.
% Vector ordering is [all x nodes at theta_1; all x nodes at theta_2; ...].
K = op.Ntheta;
nx = op.nx;
dth = op.dtheta;
u = u(:).';
rows = [];
cols = [];
vals = [];
for k = 1:K
    kp = mod(k, K) + 1;
    km = mod(k-2, K) + 1;
    a_kph = -(u(kp) - u(k))/dth;
    a_kmh = -(u(k) - u(km))/dth;
    coeff = zeros(1,K);
    coeff(k)  = coeff(k)  + max(a_kph,0)/dth - min(a_kmh,0)/dth;
    coeff(kp) = coeff(kp) + min(a_kph,0)/dth;
    coeff(km) = coeff(km) - max(a_kmh,0)/dth;
    js = find(coeff ~= 0);
    for jj = js
        rBlock = (k-1)*nx + (1:nx);
        cBlock = (jj-1)*nx + (1:nx);
        rows = [rows, rBlock]; %#ok<AGROW>
        cols = [cols, cBlock]; %#ok<AGROW>
        vals = [vals, coeff(jj)*ones(1,nx)]; %#ok<AGROW>
    end
end
T = sparse(rows, cols, vals, nx*K, nx*K);
end
