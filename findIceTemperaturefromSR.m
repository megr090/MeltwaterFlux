function [Taddiff,zetaaddiff,Br] = findIceTemperaturefromSR(strainrate,theta,H,z,Tm,Ts,K,A,n,dT,Pe)

% solve for Brinkmann number
Br = ((2*theta*H^2)./(K.*(Tm-Ts)))*(strainrate.^(n+1)./A).^(1/n);

% find critical strain rate
elatcrit = ((0.5.*Pe.^2)./(Pe-1+exp(-Pe))).^(n./(n+1)).*((K.*dT)./(A^(-1/n).*H.^2.*theta)).^(n/(n+1));

% find the thickness of the temperate zone
if strainrate > elatcrit
    zetaH = 1-(Pe/(Br))-(1/Pe).*(1+real(lambertw(-exp((-Pe.^2)./(Br)-1))));
else
    zetaH = 0;
end

zetaaddiff = zetaH.*H;

% find the temperature profile
Taddiff = zeros(size(z));
for i=1:length(z)
    if and(z(i)>=zetaH.*H,z(i)<=H)
        Taddiff(i) = Ts + dT*(Br./Pe).*(1-(z(i)./H)+(1/Pe).*exp(Pe.*(zetaH-1))-(1/Pe).*exp(Pe.*(zetaH.*H-z(i))./H));
    else
        Taddiff(i) = Tm;
    end
end

end

