%{

    Course: ELEC 4700

    PA:     10 - YEE Cell

    Date:   19 MARCH 2020

    

    Student:    Shaun Hearn

    Student ID: 100953334

%}


% The following code sets up the figure window commands
winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight



dels = 0.75;
spatialFactor = 1;


% Physical constants
c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);

% defines the time of simulation  and the frequency of solution
tSim = 200e-15;
%f = 230e12;
%f = 3e16;      % xrays
f = 4.738e14;    % red light
lambda = c_c/f;

% Sets up the spatial region of the main wave guide
xMax{1} = 20e-6;
nx{1} = 200;
ny{1} = 0.75*nx{1};


Reg.n = 1;

% Sets the permiability
mu{1} = ones(nx{1},ny{1})*c_mu_0;

% Sets the permativity
epi{1} = ones(nx{1},ny{1})*c_eps_0;

% Permitivity of the inner portion of the wave guide
% left panels
epi{1}(5:10,70:80)= c_eps_0*11.3;   
epi{1}(10:15,80:85)= c_eps_0*11.3;
epi{1}(10:15,65:70)= c_eps_0*11.3;
epi{1}(15:20,85:90)= c_eps_0*11.3;
epi{1}(15:20,60:65)= c_eps_0*11.3;
epi{1}(20:25,90:95)= c_eps_0*11.3;
epi{1}(20:25,55:60)= c_eps_0*11.3;
epi{1}(25:30,95:100)= c_eps_0*11.3;
epi{1}(25:30,50:55)= c_eps_0*11.3;
epi{1}(30:35,100:105)= c_eps_0*11.3;
epi{1}(30:35,45:50)= c_eps_0*11.3;
epi{1}(35:40,105:110)= c_eps_0*11.3;
epi{1}(35:40,40:45)= c_eps_0*11.3;
epi{1}(40:45,110:115)= c_eps_0*11.3;
epi{1}(40:45,35:40)= c_eps_0*11.3;

epi{1}(50:55,70:80)= c_eps_0*11.3;   
epi{1}(55:60,80:85)= c_eps_0*11.3;
epi{1}(55:60,65:70)= c_eps_0*11.3;
epi{1}(60:65,85:90)= c_eps_0*11.3;
epi{1}(60:65,60:65)= c_eps_0*11.3;
epi{1}(65:70,90:95)= c_eps_0*11.3;
epi{1}(65:70,55:60)= c_eps_0*11.3;
epi{1}(70:75,95:100)= c_eps_0*11.3;
epi{1}(70:75,50:55)= c_eps_0*11.3;
epi{1}(75:80,100:105)= c_eps_0*11.3;
epi{1}(75:80,45:50)= c_eps_0*11.3;
epi{1}(80:85,105:110)= c_eps_0*11.3;
epi{1}(80:85,40:45)= c_eps_0*11.3;
%epi{1}(85:90,110:115)= c_eps_0*11.3;
%epi{1}(85:90,35:40)= c_eps_0*11.3;

% Connectors
epi{1}(45:50,115:120)= c_eps_0*11.3;
epi{1}(45:50,30:35)= c_eps_0*11.3;

% Top plate
epi{1}(50:100,120:125)= c_eps_0*11.3;

% bottom plate
epi{1}(50:100,25:30)= c_eps_0*11.3;

% blinders
epi{1}(115:120,70:80)= c_eps_0*11.3;   
epi{1}(110:115,80:85)= c_eps_0*11.3;
epi{1}(110:115,65:70)= c_eps_0*11.3;
epi{1}(105:110,85:90)= c_eps_0*11.3;
epi{1}(105:110,60:65)= c_eps_0*11.3;
epi{1}(100:105,90:95)= c_eps_0*11.3;
epi{1}(100:105,55:60)= c_eps_0*11.3;
epi{1}(95:100,95:100)= c_eps_0*11.3;
epi{1}(95:100,50:55)= c_eps_0*11.3;
epi{1}(90:95,100:105)= c_eps_0*11.3;
epi{1}(90:95,45:50)= c_eps_0*11.3;
epi{1}(85:90,105:110)= c_eps_0*11.3;
epi{1}(85:90,40:45)= c_eps_0*11.3;


sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});

% Sets the grid size of the simulation and the time steps
dx = xMax{1}/nx{1};
dt = 0.25*dx/c_c;
nSteps = round(tSim/dt*2);
yMax = ny{1}*dx;
nsteps_lamda = lambda/dx;

movie = 1;
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 1.1;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];

% boundary Condiditons
bc{1}.NumS = 2; % sets the number of sources
bc{1}.s(1).xpos = nx{1}/(10) + 1; % SS position
bc{1}.s(1).type = 'ss';
bc{1}.s(1).fct = @PlaneWaveBC;

bc{1}.s(2).xpos = nx{1}/(4) + 1; % SS position
bc{1}.s(2).type = 'ss';
bc{1}.s(2).fct = @PlaneWaveBC;
% mag = -1/c_eta_0;
mag = 1;
phi = 0;
omega = f*2*pi;
betap = 0;
t0 = 30e-15;
st = 15e-15;
%st = -0.05;
s = 0;
y0 = yMax/2;
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};
bc{1}.s(2).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};

Plot.y0 = round(y0/dx);

%{
 Sets the soft source type:
    a - sets up the PML
    s - hard to H - reflects all incoming waves
    e - explicit metasurface
    m - magnetic wall (H = 0)
    p - periodic
%}
bc{1}.xm.type = 'a';
bc{1}.xp.type = 'a';
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

% Perfectly Matched Layer (PML) sizing
pml.width = 20 * spatialFactor;
pml.m = 3.5;

Reg.n  = 1; % number of regions
Reg.xoff{1} = 0;    % used for simulation surface in Yee2DEM
Reg.yoff{1} = 0;    % used for simulation surface in Yee2DEM

%{
This calls up the function RunYeeReg.  RunYeeReg initializes the E and H
fields, then calls Yee2DEM.  Yee2DEM updates the E and H fields (in
alternate time steps - E integer time steps, H half integer time steps)
and takes into consideration the Soft Source (SS) type.
%}
RunYeeReg






