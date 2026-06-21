function r = reaction_vector_1d(Kx, rho, H, extraTheta, op)
%REACTION_VECTOR_1D Build vector ordered by theta blocks.
% extraTheta is a 1-by-Ntheta vector added to H(k) in each theta block.
if nargin < 4 || isempty(extraTheta)
    extraTheta = zeros(size(H));
end
nx = op.nx;
Ktheta = op.Ntheta;
r = zeros(nx*Ktheta,1);
base = Kx(:) - rho(:);
for k = 1:Ktheta
    idx = (k-1)*nx + (1:nx);
    r(idx) = base + H(k) + extraTheta(k);
end
end
