function r = reaction_vector_1d(Kx, rho, H, extraTheta, op)
%REACTION_VECTOR_1D Build vector ordered by theta blocks.
% extraTheta is a 1-by-Ntheta vector added to H(k) in each theta block.
if nargin < 4 || isempty(extraTheta)
    extraTheta = zeros(size(H));
end
nx = op.nx;
Ktheta = op.Ntheta;
base = Kx(:) - rho(:);                    % nx x 1
ht = real(H(:).' + extraTheta(:).');       % 1 x Ktheta
rMat = bsxfun(@plus, base, ht);            % nx x Ktheta
r = reshape(rMat, nx*Ktheta, 1);
end
