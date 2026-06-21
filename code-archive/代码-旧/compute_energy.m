function E = compute_energy(Phi, V, Nx, Ny, Nz, dx, dy, dz, beta, lambda3)
   %% º∆À„Ã›∂»
    dphi_dx = zeros(Nx,Ny,Nz);
    dphi_dx(2:Nx-1,:,:)=Phi(3:Nx,:,:)-Phi(1:Nx-2,:,:)/(2*dx);
    dphi_dx(1,:,:)=Phi(2,:,:)-Phi(1,:,:)/(dx);
    dphi_dx(Nx,:,:)=Phi(Nx,:,:)-Phi(Nx-1,:,:)/(dx);
    
    dphi_dy = zeros(Nx,Ny,Nz);
    dphi_dy(2:Ny-1,:,:)=Phi(3:Ny,:,:)-Phi(1:Ny-2,:,:)/(2*dy);
    dphi_dy(1,:,:)=Phi(2,:,:)-Phi(1,:,:)/(dy);
    dphi_dy(Nx,:,:)=Phi(Ny,:,:)-Phi(Ny-1,:,:)/(dy);
    
    dphi_dz = zeros(Nx,Ny,Nz);
    dphi_dz(2:Nz-1,:,:)=Phi(3:Nz,:,:)-Phi(1:Nz-2,:,:)/(2*dz);
    dphi_dz(1,:,:)=Phi(2,:,:)-Phi(1,:,:)/(dz);
    dphi_dz(Nx,:,:)=Phi(Nz,:,:)-Phi(Nz-1,:,:)/(dz);
    
    %% 
    E_kinetic = 0.5 * (abs(dphi_dx).^2 + abs(dphi_dy).^2 + abs(dphi_dz).^2);
    E_potential = V .* abs(Phi).^2;
    E_nonlinear = 0.5 * beta * abs(Phi).^4 + (lambda3 / 4) * abs(Phi).^5;
    E = sum(E_kinetic(:) + E_potential(:) + E_nonlinear(:)) * dx * dy * dz;
end

