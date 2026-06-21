function uNew = step_phase_1d(u, H, eps, dt, op, method, lfAlpha)
%STEP_PHASE_1D Semi-implicit phase update.
Hhat = src.numerical_hamiltonian(u, op.dtheta, method, lfAlpha);
B = op.Itheta - dt * eps * op.Ltheta;
rhs = u(:) + dt * (-Hhat(:) - H(:));
uNew = (B \ rhs).';
end
