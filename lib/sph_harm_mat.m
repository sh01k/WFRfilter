%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Shoichi Koyama (koyama.shoichi@ieee.org) 2016.09.02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y] = sph_harm_mat(Nmax,theta,phi)

Y = zeros(Nmax+1,2*Nmax+1);

n = 0:Nmax;

for in = 1:(Nmax+1)
    m = -n(in):n(in);
    
    Pn = legendre(n(in),cos(theta));

    for im = 1:(2*n(in)+1)
        mp = abs(m(im));
        Pnm = Pn(mp+1);
        
        Y(in,im) = sqrt(((2.*n(in)+1)./(4.*pi)).*(factorial(n(in)-mp)./factorial(n(in)+mp))).*Pnm.*exp(1i.*m(im).*phi);
    end
end