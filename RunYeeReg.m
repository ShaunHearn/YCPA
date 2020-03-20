
global factors;
clear global MetaSurf;
clear global inputField;
global c_eps_0 c_mu_0 c_c c_eta_0
global MetaSurf;
global inputField;
global spatialFactor;
global Ez_output_t;
global Ez_output_r;

tSim = 200e-15;
k0 = 2*pi/lambda;
dt = tSim/nSteps;

X0 = xMax{1}/2;
Y0 = yMax/2;
w = xMax{1}/10;

% initialize fields to zero
Ez0{1} = zeros(nx{1},ny{1});
Hx0{1} = zeros(nx{1},ny{1}+1);
Hy0{1} = zeros(nx{1}+1,ny{1});

Ez1{1} = zeros(nx{1},ny{1});
Hx1{1} = zeros(nx{1},ny{1}+1);
Hy1{1} = zeros(nx{1}+1,ny{1});

% propagate!
[Ez_t, Ez_r, Ez,Ez1_t, Ez1_r, Ez1, Hx, Hy, Hx1, Hy1] = ...
    Yee2DEM(nx,ny,epi,mu,sigma,sigmaH,...
    xMax,tSim,nSteps,Ez0,Hx0,Hy0,Ez1,Hx1,Hy1,bc,pml,Plot,Reg,movie);

% localTitle = strcat(mainTitle, '(');
% localTitle = strcat(localTitle, 'ny=', num2str(ny{1}));
% localTitle = strcat(localTitle, ').mat');
% 
% save(localTitle, 'Ez_t', 'Ez_r', 'Ez_output_t', 'Ez_output_r', 'inputField', 'samplingLocation', 'Ez', 'Hx', 'Hy', 'dt', 'tSim', 'nSteps', 'dx', 'nx', 'ny');

