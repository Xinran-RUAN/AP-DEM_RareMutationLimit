% 构造算子 -Laplacian 对应的差分矩阵
%   周期边界条件
function[MAT_D] = Matrix_Laplacian1D_Periodic(N, dx)
e = ones(N, 1);
B = spdiags([-e, 2*e, -e], [-1, 0 ,1], N, N);
B(1, N) = - 1;
B(N, 1) = - 1;  
MAT_D = B / dx^2;
