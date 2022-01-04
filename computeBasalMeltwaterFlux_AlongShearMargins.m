function [flux,effpressure_composite,porosity_composite,T,zct] = computeBasalMeltwaterFlux_AlongShearMargins(strainrate,smb,theta,H,Ts,kappa,alpha,N0,delta,z)

n = 3;
Tm = 273; % in K
dT = Tm-Ts; % K
rho = 917; % kg m^-3
cp = 2050; % J kg^-1 K^-1
K = 2.1; % W m^-1 K^-1
A = 2.4e-24; % Pa^-3 s^-1
rho_water = 1000; % kg m^-3
rho_ice = 917; % kg m^-3;

smb = smb./(1e3*3.154e7);
Pe = (rho.*cp.*smb.*H)./(K);

[T,zct,Br] = findIceTemperaturefromSR(strainrate,theta,H,z.*H,Tm,Ts,K,A,n,dT,Pe);
zct = zct./H; % non dimensionalize

if zct > 0
    
    [porosity_outer] = findOuterPorosityProfile(-Pe,Br,zct,z,kappa,alpha);
    [effpressure_outer] = findOuterPressureProfile(Br,-Pe,porosity_outer,kappa,alpha);
    
    [porosity_inner_0] = findInnerPorosity0Profile(porosity_outer);
    [effpressure_inner] = findInnerPressureProfile(effpressure_outer,N0,alpha,kappa,-Pe,porosity_outer,z,zct,delta);
    [porosity_inner_1] = findInnerPorosity1Profile(Br,-Pe,porosity_inner_0,z,zct,delta,effpressure_outer,alpha,kappa,N0);
    [porosity_inner] = findInnerPorosityProfile(porosity_inner_0,porosity_inner_1,delta,z,zct);
    
    [effpressure_composite] = findCompositePressureProfile(effpressure_outer,effpressure_inner);
    effpressure_composite(z>zct) = 0;
    [porosity_composite] = findCompositePorosityProfile(porosity_outer,porosity_inner,Br,kappa,alpha,-Pe,zct,z);
    
    [flux] = findMeltwaterFlux(effpressure_composite,porosity_composite,alpha,kappa,z,delta);
    
else
    flux = zeros(size(T));
    effpressure_composite = zeros(size(T));
    porosity_composite = zeros(size(T));
    
end

end

