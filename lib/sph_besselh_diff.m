%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Shoichi Koyama (koyama.shoichi@ieee.org) 2016.09.02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Z] = sph_besselh_diff(N,K,X)

Z = (N.*sph_besselh(N-1,K,X)-(N+1).*sph_besselh(N+1,K,X))./(2*N+1);

end