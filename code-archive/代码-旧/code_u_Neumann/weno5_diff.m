function du = weno5_diff(u, dx)
N = length(u);
du = zeros(size(u));
eps0 = 1e-6;
u_ext = [u(end-2:end), u, u(1:3)];

for j = 1:N
    jj = j + 3;

    % ==== 正流向 WENO5 ====
    f0 = u_ext(jj-2);
    f1 = u_ext(jj-1);
    f2 = u_ext(jj);
    f3 = u_ext(jj+1);
    f4 = u_ext(jj+2);

    beta0 = (13/12)*(f0-2*f1+f2)^2 + (1/4)*(f0 -4*f1 +3*f2)^2;
    beta1 = (13/12)*(f1-2*f2+f3)^2 + (1/4)*(f1 - f3)^2;
    beta2 = (13/12)*(f2-2*f3+f4)^2 + (1/4)*(3*f2 -4*f3 +f4)^2;

    g0 = 0.1; g1 = 0.6; g2 = 0.3;
    a0 = g0 / (eps0 + beta0)^2;
    a1 = g1 / (eps0 + beta1)^2;
    a2 = g2 / (eps0 + beta2)^2;
    wsum = a0 + a1 + a2;
    w0 = a0 / wsum;
    w1 = a1 / wsum;
    w2 = a2 / wsum;

    d0 = (1/3)*f0 - (7/6)*f1 + (11/6)*f2;
    d1 = -(1/6)*f1 + (5/6)*f2 + (1/3)*f3;
    d2 = (1/3)*f2 + (5/6)*f3 - (1/6)*f4;

    du_p = (w0*d0 + w1*d1 + w2*d2)/dx;

    % ==== 负流向 WENO5 ====
    b0 = u_ext(jj+2);
    b1 = u_ext(jj+1);
    b2 = u_ext(jj);
    b3 = u_ext(jj-1);
    b4 = u_ext(jj-2);

    beta0m = (13/12)*(b0-2*b1+b2)^2 + (1/4)*(b0 -4*b1 +3*b2)^2;
    beta1m = (13/12)*(b1-2*b2+b3)^2 + (1/4)*(b1 - b3)^2;
    beta2m = (13/12)*(b2-2*b3+b4)^2 + (1/4)*(3*b2 -4*b3 +b4)^2;

    a0m = g0 / (eps0 + beta0m)^2;
    a1m = g1 / (eps0 + beta1m)^2;
    a2m = g2 / (eps0 + beta2m)^2;
    wsumm = a0m + a1m + a2m;
    w0m = a0m / wsumm;
    w1m = a1m / wsumm;
    w2m = a2m / wsumm;

    d0m = -(1/3)*b0 + (7/6)*b1 - (11/6)*b2;
    d1m = (1/6)*b1 - (5/6)*b2 - (1/3)*b3;
    d2m = -(1/3)*b2 - (5/6)*b3 + (1/6)*b4;

    du_m = (w0m*d0m + w1m*d1m + w2m*d2m)/dx;

    % ==== 对称合成 ====
    du(j) = 0.5 * (du_p + du_m);end
end

