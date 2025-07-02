% 单调格式
function flux = monotone_flux(u, dx)
    pL = WENO5_left(u, dx);
    pR = WENO5_right(u, dx);
    % grad_sq = (pL.^2) .* (pL >= 0) + (pR.^2) .* (pL < 0);
    % Osher-Sethian monotone Hamiltonian
    grad_sq = (min(pL, 0)).^2 + (max(pR, 0)).^2;
    flux = -grad_sq;
end
