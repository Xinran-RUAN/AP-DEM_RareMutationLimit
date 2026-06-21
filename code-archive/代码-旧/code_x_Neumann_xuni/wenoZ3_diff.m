function du = wenoZ3_diff(u, dx)
N = length(u);
du = zeros(size(u));
eps0 = 1e-6;
u_ext = [u(end-1:end), u, u(1:2)];

for j = 1:N
    jj = j + 2;

    %% 正向 WENO-Z3
    f0 = u_ext(jj-1);
    f1 = u_ext(jj);
    f2 = u_ext(jj+1);

    beta0 = (f1 - f0)^2;
    beta1 = (f2 - f1)^2;

    tau3 = abs(beta0 - beta1);

    d0 = 1/3;
    d1 = 2/3;

    a0 = d0 * (1 + tau3 / (eps0 + beta0));
    a1 = d1 * (1 + tau3 / (eps0 + beta1));
    wsum = a0 + a1;
    w0 = a0 / wsum;
    w1 = a1 / wsum;

    p0 = (-0.5)*f0 + (1.5)*f1;
    p1 = 0.5*f1 + 0.5*f2;

    du_p = (w0*p0 + w1*p1)/dx;

    %% 负向 WENO-Z3
    b0 = u_ext(jj+1);
    b1 = u_ext(jj);
    b2 = u_ext(jj-1);

    beta0m = (b1 - b0)^2;
    beta1m = (b2 - b1)^2;

    tau3m = abs(beta0m - beta1m);

    a0m = d0 * (1 + tau3m / (eps0 + beta0m));
    a1m = d1 * (1 + tau3m / (eps0 + beta1m));
    wsum_m = a0m + a1m;
    w0m = a0m / wsum_m;
    w1m = a1m / wsum_m;

    p0m = (-0.5)*b0 + (1.5)*b1;
    p1m = 0.5*b1 + 0.5*b2;

    du_m = -(w0m*p0m + w1m*p1m)/dx;

    %% 合成
    du(j) = du_p + du_m;
end
end

