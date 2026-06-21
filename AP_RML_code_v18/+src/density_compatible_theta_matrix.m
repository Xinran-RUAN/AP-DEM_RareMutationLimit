function B = density_compatible_theta_matrix(u, eps, op)
%DENSITY_COMPATIBLE_THETA_MATRIX Matrix for E^{-1} delta_theta^2(W E).
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
    rPlus = src.safe_exp((u(kp) - u(k))/eps);
    rMinus = src.safe_exp((u(km) - u(k))/eps);
    coeffs = [rMinus, -2, rPlus] / dth^2;
    js = [km, k, kp];
    for q = 1:3
        rBlock = (k-1)*nx + (1:nx);
        cBlock = (js(q)-1)*nx + (1:nx);
        rows = [rows, rBlock]; %#ok<AGROW>
        cols = [cols, cBlock]; %#ok<AGROW>
        vals = [vals, coeffs(q)*ones(1,nx)]; %#ok<AGROW>
    end
end
B = sparse(rows, cols, vals, nx*K, nx*K);
end
