%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Shoichi Koyama (koyama.shoichi@ieee.org) 2016.09.01
%
% 3D free-field Green's function
% USAGE "green3d(xr,yr,zr,xs,ys,zs,k)"
% INPUT
%   xr, yr, zr       Reference positions
%   xs, ys, zs       Source positions
%   k                Wave number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = green3d(xr, yr, zr, xs, ys, zs, k)
    r = sqrt((xr-xs).^2+(yr-ys).^2+(zr-zs).^2);
    G = exp(1i*k*r)./(4*pi*r);
end