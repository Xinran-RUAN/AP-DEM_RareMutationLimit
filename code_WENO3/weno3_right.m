function phix_plus = weno3_right(phi, dx)
% WENO3 ������������ع�
% phi: 1D ���飬�����ĺ���ֵ
% dx: ������
% phix_plus: �������� WENO3 ����

N = length(phi);
phix_plus = zeros(size(phi));
eps0 = 1e-1;

for i = 1:N
    % ��������
    i_m1 = mod(i-2, N) + 1; % i-1
    i_p0 = i;               % i
    i_p1 = mod(i, N) + 1;   % i+1
    i_p2 = mod(i+1, N) + 1; % i+2

    % ���
    a = (phi(i_p1) - phi(i_p0)) / dx; % ��?��_i
    b = (phi(i_p0) - phi(i_m1)) / dx; % ��?��_{i-1}
    c = (phi(i_p2) - phi(i_p1)) / dx; % ��?��_{i+1}

    % ������ģ�嵼������
    phi0 = (1/2)*a + (1/2)*b;
    phi1 = (3/2)*a - (1/2)*c;

    % �⻬��ָ��
    IS0 = (a - b)^2;
    IS1 = (a - c)^2;

    % Ȩ�أ�����Ȩ�أ�phi0 ��Ӧ 2/3��phi1 ��Ӧ 1/3��
    alpha0 = (2/3) / (eps0 + IS0)^2;
    alpha1 = (1/3) / (eps0 + IS1)^2;
    alpha_sum = alpha0 + alpha1;
    w0 = alpha0 / alpha_sum;
    w1 = alpha1 / alpha_sum;

    % �ϳɵ���
    phix_plus(i) = w0 * phi0 + w1 * phi1;
end
end
