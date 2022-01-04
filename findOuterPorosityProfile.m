function [porosity_outer] = findOuterPorosityProfile(Pe,Br,zct,z,kappa,alpha)

options = optimoptions('fsolve','Diagnostics','off','Display','off','MaxIterations',500);
fun = @(x)PorosityOuterEqn(x,Pe,Br,zct,z,kappa,alpha);
x0 = 1.*ones(size(z))+z;
[x] = fsolve(fun,x0,options);
porosity_outer = real(x);

porosity_outer(z>zct) = 0;

end