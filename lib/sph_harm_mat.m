%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Shoichi Koyama 2015.07.31
% Generate matrix of spherical harmonic functions
% USAGE "sph_harm_mat(Nmax,theta,phi)"
% INPUT
%   Nmax     Maximum order  
%   theta    Zenith angle
%   phi      Azimuth angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y] = sph_harm_mat(Nmax,theta,phi)

Y = zeros(Nmax+1,2*Nmax+1);

n = 0:1:Nmax;
size_n = length(n);

for in=1:size_n
    m = -n(in):1:n(in);
    size_m = length(m);
    
    %Associated Legendre function
    Pnm = legendre(n(in),cos(theta));
    Pnm_vec = [Pnm(n(in)+1:-1:2); Pnm].';
    
    sgn = [(-1).^abs(m(1:n(in))),ones(1,n(in)+1)];
    Y(in,1:size_m) = sgn.*sqrt((2*n(in)+1)/(4*pi).*factorial(n(in)-abs(m))./factorial(n(in)+abs(m))).*Pnm_vec.*exp(1i*m*phi);
end

end