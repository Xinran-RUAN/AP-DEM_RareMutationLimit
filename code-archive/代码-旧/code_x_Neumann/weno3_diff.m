function du = weno3_diff(u, dx)
N = length(u);
du = zeros(size(u));
eps0 = 1e-6;
u_ext = [u(end-1:end), u, u(1:2)];

for j = 1:N
    jj = j + 2;

    % ==== 正流向 ====
    f_im1 = u_ext(jj-1);
    f_i   = u_ext(jj);
    f_ip1 = u_ext(jj+1);

    beta0 = (f_i - f_im1)^2;
    beta1 = (f_ip1 - f_i)^2;

    d0 = 1/3;
    d1 = 2/3;

    a0 = d0 / (eps0 + beta0)^2;
    a1 = d1 / (eps0 + beta1)^2;
    wsum = a0 + a1;
    w0 = a0 / wsum;
    w1 = a1 / wsum;

    p0 = (-f_im1 + 3*f_i)/2;
    p1 = (f_i + f_ip1)/2;

    du_p = (w0*p0 + w1*p1)/dx;

    % ==== 负流向 ====
    beta0m = (f_i - f_ip1)^2;
    beta1m = (f_im1 - f_i)^2;

    a0m = d0 / (eps0 + beta0m)^2;
    a1m = d1 / (eps0 + beta1m)^2;
    wsum_m = a0m + a1m;
    w0m = a0m / wsum_m;
    w1m = a1m / wsum_m;

    p0m = (-f_ip1 + 3*f_i)/2;
    p1m = (f_i + f_im1)/2;

    du_m = -(w0m*p0m + w1m*p1m)/dx;

    % 合成
    du(j) = du_p + du_m;
end
end

