function phix_plus = weno5_right(phi, dx)
% WENO5 ������������ع�
% phi: 1D ���飬�����ĺ���ֵ
% dx: ������
% phix_plus: �������� WENO5 ����

N = length(phi);
phix_plus = zeros(size(phi));
eps0 = 1e-6;

for i = 1:N
    % ��������
    i_m2 = mod(i-3, N) + 1; % i-2
    i_m1 = mod(i-2, N) + 1; % i-1
    i_p1 = mod(i, N) + 1;   % i+1
    i_p2 = mod(i+1, N) + 1; % i+2
    i_p3 = mod(i+2, N) + 1; % i+3

    % ���
    a = (phi(i_p2) - phi(i_p1)) / dx; % ��?��_{i+1}
    b = (phi(i_p1) - phi(i)) / dx;    % ��?��_i
    c = (phi(i) - phi(i_m1)) / dx;    % ��?��_{i-1}
    d = (phi(i_m1) - phi(i_m2)) / dx; % ��?��_{i-2}

    % ������ģ�嵼��
    phi0 = (1/3)*( (phi(i_p3)-phi(i_p2))/dx ) - (7/6)*a + (11/6)*b;
    phi1 = -(1/6)*a + (5/6)*b + (1/3)*c;
    phi2 = (1/3)*b + (5/6)*c - (1/6)*d;

    % �⻬��ָ��
    IS0 = 13*((phi(i_p3)-phi(i_p2))/dx - a)^2 + 3*((phi(i_p3)-phi(i_p2))/dx - 3*a)^2;
    IS1 = 13*(a - b)^2 + 3*(a - b)^2;
    IS2 = 13*(b - c)^2 + 3*(3*b - c)^2;

    % Ȩ��
    alpha0 = 1 / (eps0 + IS0)^2;
    alpha1 = 6 / (eps0 + IS1)^2;
    alpha2 = 3 / (eps0 + IS2)^2;
    alpha_sum = alpha0 + alpha1 + alpha2;
    w0 = alpha0 / alpha_sum;
    w1 = alpha1 / alpha_sum;
    w2 = alpha2 / alpha_sum;

    % �ϳɵ���
    phix_plus(i) = w0*phi0 + w1*phi1 + w2*phi2;
end
end
