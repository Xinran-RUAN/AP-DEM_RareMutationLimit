function du = wenoZ5_diff(u, dx)
N = length(u);
du = zeros(size(u));
eps0 = 1e-6; 
u_ext = [u(end-2:end), u, u(1:3)];

for j = 1:N
    jj = j + 3;

    %% 正向 WENO-Z
    f0 = u_ext(jj-2);
    f1 = u_ext(jj-1);
    f2 = u_ext(jj);
    f3 = u_ext(jj+1);
    f4 = u_ext(jj+2);

    beta0 = (13/12)*(f0 - 2*f1 + f2)^2 + (1/4)*(f0 - 4*f1 + 3*f2)^2;
    beta1 = (13/12)*(f1 - 2*f2 + f3)^2 + (1/4)*(f1 - f3)^2;
    beta2 = (13/12)*(f2 - 2*f3 + f4)^2 + (1/4)*(3*f2 - 4*f3 + f4)^2;

    tau5 = abs(beta0 - beta2);

    d0 = 1/10; d1 = 6/10; d2 = 3/10;
    a0 = d0 * (1 + tau5 / (eps0 + beta0));
    a1 = d1 * (1 + tau5 / (eps0 + beta1));
    a2 = d2 * (1 + tau5 / (eps0 + beta2));
    wsum = a0 + a1 + a2;
    w0 = a0 / wsum;
    w1 = a1 / wsum;
    w2 = a2 / wsum;

    d0f = (1/3)*f0 - (7/6)*f1 + (11/6)*f2;
    d1f = -(1/6)*f1 + (5/6)*f2 + (1/3)*f3;
    d2f = (1/3)*f2 + (5/6)*f3 - (1/6)*f4;

    du_p = (w0*d0f + w1*d1f + w2*d2f)/dx;

    %% 负向 WENO-Z
    b0 = u_ext(jj+2);
    b1 = u_ext(jj+1);
    b2 = u_ext(jj);
    b3 = u_ext(jj-1);
    b4 = u_ext(jj-2);

    beta0m = (13/12)*(b0 - 2*b1 + b2)^2 + (1/4)*(b0 - 4*b1 + 3*b2)^2;
    beta1m = (13/12)*(b1 - 2*b2 + b3)^2 + (1/4)*(b1 - b3)^2;
    beta2m = (13/12)*(b2 - 2*b3 + b4)^2 + (1/4)*(3*b2 - 4*b3 + b4)^2;

    tau5m = abs(beta0m - beta2m);

    a0m = d0 * (1 + tau5m / (eps0 + beta0m));
    a1m = d1 * (1 + tau5m / (eps0 + beta1m));
    a2m = d2 * (1 + tau5m / (eps0 + beta2m));
    wsumm = a0m + a1m + a2m;
    w0m = a0m / wsumm;
    w1m = a1m / wsumm;
    w2m = a2m / wsumm;

    d0fm = -(1/3)*b0 + (7/6)*b1 - (11/6)*b2;
    d1fm = (1/6)*b1 - (5/6)*b2 - (1/3)*b3;
    d2fm = -(1/3)*b2 - (5/6)*b3 + (1/6)*b4;

    du_m = (w0m*d0fm + w1m*d1fm + w2m*d2fm)/dx;

    %% 合成
    du(j) = du_p + du_m;
end
end

