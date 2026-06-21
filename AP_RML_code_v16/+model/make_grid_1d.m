function grid = make_grid_1d(par)
%MAKE_GRID_1D Endpoint grid in x and periodic cell grid in theta.
grid.Nx = par.Nx;
grid.Ntheta = par.Ntheta;
grid.x = linspace(0, 1, par.Nx + 1).';
grid.dx = 1 / par.Nx;
grid.theta = (0:par.Ntheta-1) / par.Ntheta;
grid.dtheta = 1 / par.Ntheta;
end
