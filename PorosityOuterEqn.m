function [func] = PorosityOuterEqn(phi,Pe,Br,zct,z,kappa,alpha)

func = Pe.*phi-kappa.*phi.^alpha-Br.*(z-zct);

end