function [effpressure_inner] = findInnerPressureProfile(effpressure_outer,N0,alpha,kappa,Pe,porosity_inner_0,z,zct,delta)


lambda = sqrt(((Pe-alpha.*kappa.*porosity_inner_0.^(alpha-1))).*porosity_inner_0./(Pe.*kappa.*porosity_inner_0.^alpha));
lambda = real(lambda);

effpressure_inner = effpressure_outer(1)+(N0-effpressure_outer(1)).*exp(-lambda.*z.*(1./(delta.^(1/2))));
effpressure_inner(z>zct) = 0;

end

