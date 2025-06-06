function some_try
% --------- �������� ----------
tic 
N = 10000;                   % ���������ɸ�Ϊ����
h = 1 / (N + 1);            % ���񲽳�
x = linspace(0, 1, N+2);    % �����߽��
e = ones(N,1);

% --------- �������޲�־���Dirichlet �߽磩 ----------
A = spdiags([-e 2*e -e], -1:1, N, N) / h^2;

% --------- Shift-and-Invert ���ݷ�����С����ֵ ----------
sigma = 0;                          % shift��Ŀ����Сֵ������
S = A - sigma * speye(N);           % shifted ����
[L,U] = lu(S);                      % LU �ֽ⣨һ���ԣ�

w = rand(N,1);                      % ��ʼ����
w = w / norm(w);
max_iter = 100;
tol = 1e-10;

for k = 1:max_iter
    w_old = w;
    w = U \ (L \ w);               % �� (A - sigma*I)x = b
    w = w / norm(w);               % ��λ��
    if norm(w - w_old) < tol
        break;
    end
end

lambda = (w' * A * w) / (w' * w);  % Rayleigh �̹�������ֵ

% --------- ������ ----------
fprintf("��С����ֵ����: lambda = %.12f\n", lambda);
fprintf("��������: %d\n", k);
toc
end

