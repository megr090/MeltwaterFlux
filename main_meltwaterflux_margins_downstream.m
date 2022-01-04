%% Computing meltwater flux from shear margins

% Here we estimate effective pressure and effective porosity, and then
% meltwater flux from shear margins in Antarctica

clear all;
close all;

% choose ice stream
icestream = 'byrd'; % 'pig', 'byrd', 'bindschadler'

% top-level parameters
kappa = 0.52;
alpha = 2.33;
delta = 0.001;
num_pts = 100;
K = 2.1; % W m^-1 K^-1
rho_w = 1000; % kg m^-3
rho_ice = 917; % kg m^-3;
epsilon = 0.01;
L = 3.34e5;
n = 3;

[SRvec,SMBvec,Hvec,Tsvec,Vvec] = readData_AlongShearMargin(icestream);

% find basal boundary condition for N
N0_dim = 20e3;
dT = 273-mean(Tsvec);
refviscosity = 0.5.*2.4e-24.^(-1/n).*mean(SRvec).^((1/n)-1);
refN_dimcomp = refviscosity.*K.*dT./epsilon./rho_w./L./mean(Hvec).^2;
N0 = N0_dim./refN_dimcomp;

x = ones(length(SRvec),1);
x(end) = 0;
for i=length(SRvec)-1:-1:1
    x(i) = x(i+1)+500;
end

% initialize fields
flux = zeros(num_pts,length(SRvec));
porosity_composite = zeros(num_pts,length(SRvec));
pressure_composite = zeros(num_pts,length(SRvec));
T = zeros(num_pts,length(SRvec));
A = zeros(num_pts,length(SRvec));
viscosity = zeros(num_pts,length(SRvec));
N_dimcomp = zeros(num_pts,length(SRvec));
N_dim = zeros(num_pts,length(SRvec));
flux_dim = zeros(num_pts,length(SRvec));
zct = zeros(length(SRvec),1);

for i=1:length(SRvec)
    
    H = Hvec(i);
    z = linspace(0,H,num_pts);
    z = z./H; % non dimensionalize
    
    Ts = Tsvec(i);
    
    strainrate = SRvec(i);
    smb = SMBvec(i);
    
    theta = 1;
    
    [flux(:,i),pressure_composite(:,i),porosity_composite(:,i),T(:,i),zct(i)] = computeBasalMeltwaterFlux_AlongShearMargins(strainrate,smb,theta,H,Ts,kappa,alpha,N0,delta,z);
    dT = 273-mean(Tsvec);
    viscosity(:,i) = 0.5.*2.4e-24.^(-1/n).*mean(SRvec).^((1/n)-1);
    N_dimcomp(:,i) = viscosity(:,i).*K.*dT./epsilon./rho_w./L./mean(Hvec).^2;
    N_dim(:,i) = N_dimcomp(:,i).*pressure_composite(:,i);
    
    [flux_dim(:,i)] = (K.*dT./(rho_w.*L.*mean(Hvec))).*flux(:,i);
    flux_dimcomp = (K.*dT./(rho_w.*L.*mean(Hvec)));
    
    fprintf('Iteration %d of %d done \n',i,length(SRvec));
end

%% Plotting

[done] = plotMeltwaterFluxAlongMargins(x,flux,pressure_composite,porosity_composite,N_dim,flux_dim,N_dimcomp,SRvec,Hvec,SMBvec,Tsvec,epsilon,N0,zct,flux_dimcomp);
