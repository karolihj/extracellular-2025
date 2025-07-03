function B = set_up_sub_matrix(G, mesh, M)
%B = set_up_sub_matrix(G, mesh, M)

% Load parameters
N = G.N;
Nx = G.Nx;
Ny = G.Ny;

% Load grid point names
e_lsw = mesh.e_lsw;
e_lse = mesh.e_lse;
e_lnw = mesh.e_lnw;
e_lne = mesh.e_lne;
e_hsw = mesh.e_hsw;
e_hse = mesh.e_hse;
e_hnw = mesh.e_hnw;
e_hne = mesh.e_hne;

e_hw = mesh.e_hw;
e_he = mesh.e_he;
e_hs = mesh.e_hs;
e_hn = mesh.e_hn;
e_lw = mesh.e_lw;
e_le = mesh.e_le;
e_ls = mesh.e_ls;
e_ln = mesh.e_ln;
e_ne = mesh.e_ne;
e_sw = mesh.e_sw;
e_se = mesh.e_se;
e_nw = mesh.e_nw;

e_w = mesh.e_w;
e_e = mesh.e_e;
e_s = mesh.e_s;
e_n = mesh.e_n;
e_h = mesh.e_h;
e_l = mesh.e_l;

e = mesh.e;

vec = zeros(N, 1);
vec_kp = zeros(N, 1);
vec_km = zeros(N, 1);
vec_jp = zeros(N, 1);
vec_jm = zeros(N, 1);
vec_qp = zeros(N, 1);
vec_qm = zeros(N, 1);


%%%%%%% BOUNDARY %%%%%%%

% 1a) Set up rows for the low boundary
index = e_l;
vec(index) = M.px(index) + M.mx(index) + M.py(index) + M.my(index) + ...
    M.pz(index);
vec_kp(index+1) = -M.px(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_jm(index-Nx) = -M.my(index);
vec_qp(index+Nx*Ny) = -M.pz(index);

% 2a) Set up rows for the high boundary
index = e_h;
vec(index) = M.px(index) + M.mx(index) + M.py(index) + M.my(index) + ...
    M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_jm(index-Nx) = -M.my(index);
vec_qm(index-Nx*Ny) = -M.mz(index);

% 3a) Set up rows for the south boundary
index = e_s;
vec(index) = M.px(index) + M.mx(index) + M.py(index) + ...
    M.pz(index) + M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_qp(index+Nx*Ny) = -M.pz(index);
vec_qm(index-Nx*Ny) = -M.mz(index);


% 4a) Set up rows for the north boundary
index = e_n;
vec(index) = M.px(index) + M.mx(index) + M.my(index) + ...
    M.pz(index) + M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_km(index-1) = -M.mx(index);
vec_jm(index-Nx) = -M.my(index);
vec_qp(index+Nx*Ny) = -M.pz(index);
vec_qm(index-Nx*Ny) = -M.mz(index);

% 5a) Set up rows for the left boundary
index = e_w;
vec(index) = M.px(index) + M.py(index) + M.my(index) + ...
    M.pz(index) + M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_jp(index+Nx) = -M.py(index);
vec_jm(index-Nx) = -M.my(index);
vec_qp(index+Nx*Ny) = -M.pz(index);
vec_qm(index-Nx*Ny) = -M.mz(index);

% 6a) Set up rows for the right boundary
index = e_e;
vec(index) = M.mx(index) + M.py(index) + M.my(index) + ...
    M.pz(index) + M.mz(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_jm(index-Nx) = -M.my(index);
vec_qp(index+Nx*Ny) = -M.pz(index);
vec_qm(index-Nx*Ny) = -M.mz(index);

% 1b) Set up rows for the high left boundary
index = e_hw;
vec(index) = M.px(index) + M.py(index) + M.my(index) + ...
    M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_jp(index+Nx) = -M.py(index);
vec_jm(index-Nx) = -M.my(index);
vec_qm(index-Nx*Ny) = -M.mz(index);


% 2b) Set up rows for the high right boundary
index = e_he;
vec(index) = M.px(index) + M.py(index) + M.my(index) + ...
    M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_jp(index+Nx) = -M.py(index);
vec_jm(index-Nx) = -M.my(index);
vec_qm(index-Nx*Ny) = -M.mz(index);

% 3b) Set up rows for the high south boundary
index = e_hs;
vec(index) = M.px(index) + M.mx(index) + M.py(index) + M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_qm(index-Nx*Ny) = -M.mz(index);


% 4b) Set up rows for the high north boundary
index = e_hn;
vec(index) = M.px(index) + M.mx(index) + M.my(index) + ...
    M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_km(index-1) = -M.mx(index);
vec_jm(index-Nx) = -M.my(index);
vec_qm(index-Nx*Ny) = -M.mz(index);


% 5b) Set up rows for the low left boundary
index = e_lw;
vec(index) = M.px(index) + M.py(index) + M.my(index) + ...
    M.pz(index);
vec_kp(index+1) = -M.px(index);
vec_jp(index+Nx) = -M.py(index);
vec_jm(index-Nx) = -M.my(index);
vec_qp(index+Nx*Ny) = -M.pz(index);


% 6b) Set up rows for the low right boundary
index = e_le;
vec(index) = M.mx(index) + M.py(index) + M.my(index) + ...
    M.pz(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_jm(index-Nx) = -M.my(index);
vec_qp(index+Nx*Ny) = -M.pz(index);


% 7b) Set up rows for the low south boundary
index = e_ls;
vec(index) = M.px(index) + M.mx(index) + M.py(index) + M.pz(index);
vec_kp(index+1) = -M.px(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_qp(index+Nx*Ny) = -M.pz(index);


% 8b) Set up rows for the low north boundary
index = e_ln;
vec(index) = M.px(index) + M.mx(index) + M.my(index) + ...
    M.pz(index);
vec_kp(index+1) = -M.px(index);
vec_km(index-1) = -M.mx(index);
vec_jm(index-Nx) = -M.my(index);
vec_qp(index+Nx*Ny) = -M.pz(index);


% 9b) Set up rows for the north left boundary
index = e_nw;
vec(index) = M.px(index) + M.my(index) + ...
    M.pz(index) + M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_jm(index-Nx) = -M.my(index);
vec_qp(index+Nx*Ny) = -M.pz(index);
vec_qm(index-Nx*Ny) = -M.mz(index);


% 10b) Set up rows for the north right boundary
index = e_ne;
vec(index) = M.mx(index) + M.my(index) + ...
    M.pz(index) + M.mz(index);
vec_km(index-1) = -M.mx(index);
vec_jm(index-Nx) = -M.my(index);
vec_qp(index+Nx*Ny) = -M.pz(index);
vec_qm(index-Nx*Ny) = -M.mz(index);


% 11b) Set up rows for the south left boundary
index = e_sw;
vec(index) = M.px(index) + M.py(index) + ...
    M.pz(index) + M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_jp(index+Nx) = -M.py(index);
vec_qp(index+Nx*Ny) = -M.pz(index);
vec_qm(index-Nx*Ny) = -M.mz(index);


% 12b) Set up rows for the south east boundary
index = e_se;
vec(index) = M.mx(index) + M.py(index) + ...
    M.pz(index) + M.mz(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_qp(index+Nx*Ny) = -M.pz(index);
vec_qm(index-Nx*Ny) = -M.mz(index);


% 1c) Set up rows for the lower, south, left boundary
index = e_lsw;
vec(index) = M.px(index) + M.py(index) + M.pz(index);
vec_kp(index+1) = -M.px(index);
vec_jp(index+Nx) = -M.py(index);
vec_qp(index+Nx*Ny) = -M.pz(index);


% 2c) Set up rows for the lower, south, east boundary
index = e_lse;
vec(index) = M.mx(index) + M.py(index) + M.pz(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_qp(index+Nx*Ny) = -M.pz(index);


% 3c) Set up rows for the lower, north, left boundary
index = e_lnw;
vec(index) = M.px(index) + M.my(index) + M.pz(index);
vec_kp(index+1) = -M.px(index);
vec_jm(index-Nx) = -M.my(index);
vec_qp(index+Nx*Ny) = -M.pz(index);

% 4c) Set up rows for the lower, north, right boundary
index = e_lne;
vec(index) = M.mx(index) + M.my(index) + M.pz(index);
vec_km(index-1) = -M.mx(index);
vec_jm(index-Nx) = -M.my(index);
vec_qp(index+Nx*Ny) = -M.pz(index);


% 5c) Set up rows for the higher, south, left boundary
index = e_hsw;
vec(index) = M.px(index) + M.py(index) + M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_jp(index+Nx) = -M.py(index);
vec_qm(index-Nx*Ny) = -M.mz(index);

% 6c) Set up rows for the higher, south, east boundary
index = e_hse;
vec(index) = M.mx(index) + M.py(index) + M.mz(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_qm(index-Nx*Ny) = -M.mz(index);


% 7c) Set up rows for the higher, north, left boundary
index = e_hnw;
vec(index) = M.px(index) + M.my(index) + M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_jm(index-Nx) = -M.my(index);
vec_qm(index-Nx*Ny) = -M.mz(index);


% 8c) Set up rows for the higher, north, right boundary
index = e_hne;
vec(index) = M.mx(index) + M.my(index) + M.mz(index);
vec_km(index-1) = -M.mx(index);
vec_jm(index-Nx) = -M.my(index);
vec_qm(index-Nx*Ny) = -M.mz(index);

  

%%%%%%% INNER DOMAIN %%%%%%%
index = e;
vec(index) = M.px(index) + M.mx(index) + M.py(index) + M.my(index) + ...
    M.pz(index) + M.mz(index);
vec_kp(index+1) = -M.px(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_jm(index-Nx) = -M.my(index);
vec_qp(index+Nx*Ny) = -M.pz(index);
vec_qm(index-Nx*Ny) = -M.mz(index);

%%%%%%% SET UP THE MATRIX %%%%%%%
B = spdiags(vec, 0, N, N);
B = B + spdiags(vec_kp, 1, N, N);
B = B + spdiags(vec_km, -1, N, N);
B = B + spdiags(vec_jp, Nx, N, N);
B = B + spdiags(vec_jm, -Nx, N, N);
B = B + spdiags(vec_qp, Nx*Ny, N, N);
B = B + spdiags(vec_qm, -Nx*Ny, N, N);

end

