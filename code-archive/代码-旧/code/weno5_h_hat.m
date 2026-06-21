function grad_sq = weno5_h_hat(u, dx)

pL = weno5_left(u, dx);
pR = weno5_right(u, dx);

grad_sq = RF_flux(pR, pL);

end

