function [porosity_inner_1] = findInnerPorosity1Profile(Br,Pe,porosity_inner_0,z,zct,delta,effpressure_outer,alpha,kappa,N0)


lambda = sqrt(((Pe-alpha.*kappa.*porosity_inner_0.^(alpha-1))).*porosity_inner_0./(Pe.*kappa.*porosity_inner_0.^alpha));
lambda = real(lambda);
porosity_inner_1 = z.*((Br-porosity_inner_0.*effpressure_outer(1))./(delta.^0.5.*Pe))+...
    (porosity_inner_0./(lambda.*Pe)).*(N0-effpressure_outer(1)).*exp(-lambda.*(z)./delta.^0.5);
porosity_inner_1(z>zct) = 0;

end