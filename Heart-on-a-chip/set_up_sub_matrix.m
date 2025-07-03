function B = set_up_sub_matrix(G, mesh, M)
%B = set_up_sub_matrix(G, mesh, M)

% Load parameters
N = G.N;
Nx = G.Nx;

% Load grid point names
e_ne = mesh.e_ne;
e_sw = mesh.e_sw;
e_se = mesh.e_se;
e_nw = mesh.e_nw;

e_w = mesh.e_w;
e_e = mesh.e_e;
e_s = mesh.e_s;
e_n = mesh.e_n;

e = mesh.e;

vec = zeros(N, 1);
vec_kp = zeros(N, 1);
vec_km = zeros(N, 1);
vec_jp = zeros(N, 1);
vec_jm = zeros(N, 1);


%%%%%%% BOUNDARY %%%%%%%

% 1a) Set up rows for the south boundary
index = e_s;
vec(index) = M.px(index) + M.mx(index) + M.py(index);
vec_kp(index+1) = -M.px(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);


% 2a) Set up rows for the north boundary
index = e_n;
vec(index) = M.px(index) + M.mx(index) + M.my(index);
vec_kp(index+1) = -M.px(index);
vec_km(index-1) = -M.mx(index);
vec_jm(index-Nx) = -M.my(index);

% 3a) Set up rows for the left boundary
index = e_w;
vec(index) = M.px(index) + M.py(index) + M.my(index);
vec_kp(index+1) = -M.px(index);
vec_jp(index+Nx) = -M.py(index);
vec_jm(index-Nx) = -M.my(index);

% 4a) Set up rows for the right boundary
index = e_e;
vec(index) = M.mx(index) + M.py(index) + M.my(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_jm(index-Nx) = -M.my(index);

% 1b) Set up rows for the north left boundary
index = e_nw;
vec(index) = M.px(index) + M.my(index);
vec_kp(index+1) = -M.px(index);
vec_jm(index-Nx) = -M.my(index);


% 2b) Set up rows for the north right boundary
index = e_ne;
vec(index) = M.mx(index) + M.my(index);
vec_km(index-1) = -M.mx(index);
vec_jm(index-Nx) = -M.my(index);

% 3b) Set up rows for the south left boundary
index = e_sw;
vec(index) = M.px(index) + M.py(index);
vec_kp(index+1) = -M.px(index);
vec_jp(index+Nx) = -M.py(index);

% 4b) Set up rows for the south east boundary
index = e_se;
vec(index) = M.mx(index) + M.py(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);

%%%%%%% INNER DOMAIN %%%%%%%
index = e;
vec(index) = M.px(index) + M.mx(index) + M.py(index) + M.my(index);
vec_kp(index+1) = -M.px(index);
vec_km(index-1) = -M.mx(index);
vec_jp(index+Nx) = -M.py(index);
vec_jm(index-Nx) = -M.my(index);

%%%%%%% SET UP THE MATRIX %%%%%%%
B = spdiags(vec, 0, N, N);
B = B + spdiags(vec_kp, 1, N, N);
B = B + spdiags(vec_km, -1, N, N);
B = B + spdiags(vec_jp, Nx, N, N);
B = B + spdiags(vec_jm, -Nx, N, N);

end

