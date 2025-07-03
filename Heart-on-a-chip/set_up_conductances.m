function [Mm, Me] = set_up_conductances(G)
%[Mm, Me] = set_up_conductances(G)


delta_mx = (G.delta_m(2:end,:)+G.delta_m(1:end-1,:))/2;
delta_ex = (G.delta_e(2:end,:)+G.delta_e(1:end-1,:))/2;
delta_my = (G.delta_m(:,2:end)+G.delta_m(:,1:end-1))/2;
delta_ey = (G.delta_e(:,2:end)+G.delta_e(:,1:end-1))/2;

MMx = (1./(G.lx./(delta_mx*G.ly*G.lz*G.sigma_m) + G.Rg_mx));
MMy = (1./(G.ly./(delta_my*G.lx*G.lz*G.sigma_m) + G.Rg_my));
MEx = delta_ex*(G.ly*G.lz*G.sigma_e)/G.lx;
MEy = delta_ey*(G.lx*G.lz*G.sigma_e)/G.ly;


Mm.mx = reshape([zeros(1,G.Ny); MMx], G.N, 1);
Mm.px = reshape([MMx; zeros(1,G.Ny)], G.N, 1);
Me.mx = reshape([zeros(1,G.Ny); MEx], G.N, 1);
Me.px = reshape([MEx; zeros(1,G.Ny)], G.N, 1);

Mm.my = reshape([zeros(G.Nx, 1), MMy], G.N, 1);
Mm.py = reshape([MMy, zeros(G.Nx, 1)], G.N, 1);
Me.my = reshape([zeros(G.Nx, 1), MEy], G.N, 1);
Me.py = reshape([MEy, zeros(G.Nx, 1)], G.N, 1);

end