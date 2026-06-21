function pR = weno5_right_derivative(u, dx)
%WENO5_RIGHT_DERIVATIVE High-order right-biased point derivative reconstruction.
u = u(:);
eps0 = 1e-6;
N = length(u);
ue = [u(end-2:end); u; u(1:3)];
pR = zeros(N, 1);
for j = 1:N
    jj = j + 3;
    vppp = ue(jj+3); vpp = ue(jj+2); vp = ue(jj+1);
    v = ue(jj); vm = ue(jj-1); vmm = ue(jj-2);
    a = vppp - 2*vpp + vp;
    b = vpp - 2*vp + v;
    c = vp - 2*v + vm;
    d = v - 2*vm + vmm;
    IS0 = 13*(a-b)^2 + 3*(a-3*b)^2;
    IS1 = 13*(b-c)^2 + 3*(b+c)^2;
    IS2 = 13*(c-d)^2 + 3*(3*c-d)^2;
    a0 = 1/(eps0 + IS0)^2;
    a1 = 6/(eps0 + IS1)^2;
    a2 = 3/(eps0 + IS2)^2;
    wsum = a0 + a1 + a2;
    w0 = a0/wsum;
    w2 = a2/wsum;
    pR(j) = (-(vm-vmm) + 7*(v-vm) + 7*(vp-v) - (vpp-vp))/12 ...
            + w0*(a-2*b+c)/3 + (w2-1/2)*(b-2*c+d)/6;
end
pR = (pR/dx).';
end
