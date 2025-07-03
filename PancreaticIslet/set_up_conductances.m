function [Mm, Me] = set_up_conductances(G)
%[Mm, Me] = set_up_conductances(G)


delta_mx = (G.delta_m(2:end,:,:)+G.delta_m(1:end-1,:,:))/2;
delta_ex = (G.delta_e(2:end,:,:)+G.delta_e(1:end-1,:,:))/2;
delta_my = (G.delta_m(:,2:end,:)+G.delta_m(:,1:end-1,:))/2;
delta_ey = (G.delta_e(:,2:end,:)+G.delta_e(:,1:end-1,:))/2;
delta_mz = (G.delta_m(:,:,2:end)+G.delta_m(:,:,1:end-1))/2;
delta_ez = (G.delta_e(:,:,2:end)+G.delta_e(:,:,1:end-1))/2;

MMx = (1./(G.lx./(delta_mx*G.ly*G.lz*G.sigma_m) + G.Rg_mx));
MMy = (1./(G.ly./(delta_my*G.lx*G.lz*G.sigma_m) + G.Rg_my));
MMz = (1./(G.lz./(delta_mz*G.lx*G.ly*G.sigma_m) + G.Rg_mz));
MEx = delta_ex*(G.ly*G.lz*G.sigma_e)/G.lx;
MEy = delta_ey*(G.lx*G.lz*G.sigma_e)/G.ly;
MEz = delta_ez*(G.lx*G.ly*G.sigma_e)/G.lz;

Mm.mx = zeros(G.Nx, G.Ny, G.Nz);
Mm.mx(2:end,:,:) = MMx;
Mm.mx = reshape(Mm.mx, G.N, 1);
Mm.px = zeros(G.Nx, G.Ny, G.Nz);
Mm.px(1:end-1,:,:) = MMx;
Mm.px = reshape(Mm.px, G.N, 1);
Me.mx = zeros(G.Nx, G.Ny, G.Nz);
Me.mx(2:end,:,:) = MEx;
Me.mx = reshape(Me.mx, G.N, 1);
Me.px = zeros(G.Nx, G.Ny, G.Nz);
Me.px(1:end-1,:,:) = MEx;
Me.px = reshape(Me.px, G.N, 1);

Mm.my = zeros(G.Nx, G.Ny, G.Nz);
Mm.my(:,2:end,:) = MMy;
Mm.my = reshape(Mm.my, G.N, 1);
Mm.py = zeros(G.Nx, G.Ny, G.Nz);
Mm.py(:,1:end-1,:) = MMy;
Mm.py = reshape(Mm.py, G.N, 1);
Me.my = zeros(G.Nx, G.Ny, G.Nz);
Me.my(:,2:end,:) = MEy;
Me.my = reshape(Me.my, G.N, 1);
Me.py = zeros(G.Nx, G.Ny, G.Nz);
Me.py(:,1:end-1,:) = MEy;
Me.py = reshape(Me.py, G.N, 1);

Mm.mz = zeros(G.Nx, G.Ny, G.Nz);
Mm.mz(:,:,2:end) = MMz;
Mm.mz = reshape(Mm.mz, G.N, 1);
Mm.pz = zeros(G.Nx, G.Ny, G.Nz);
Mm.pz(:,:,1:end-1) = MMz;
Mm.pz = reshape(Mm.pz, G.N, 1);
Me.mz = zeros(G.Nx, G.Ny, G.Nz);
Me.mz(:,:,2:end) = MEz;
Me.mz = reshape(Me.mz, G.N, 1);
Me.pz = zeros(G.Nx, G.Ny, G.Nz);
Me.pz(:,:,1:end-1) = MEz;
Me.pz = reshape(Me.pz, G.N, 1);



end