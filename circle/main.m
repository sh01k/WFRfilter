%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Shoichi Koyama (koyama.shoichi@ieee.org) 2016.09.02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables;
close all;

addpath('../lib');

% set(0,'defaultAxesFontSize',18);
% set(0,'defaultAxesFontName','Times');
% set(0,'defaultTextFontSize',18);
% set(0,'defaultTextFontName','Times');

%% Parameters

WindowLocation;

%Sound speed [m/s]
c=340.29;

%% Parameters: Reproduced area
%Size of reproduced area [m]
lenX=3.3;
lenY=3.3;

%Interval of reproduced area [m]
dx=0.015;
dy=0.015;

%Number of samples
Nx=round(lenX/dx);
Ny=round(lenY/dy);

%Pisitions
xr = ((0:Nx-1)-Nx/2).*dx;
yr = ((0:Ny-1)-Ny/2).*dy;
zr = 0.0;

%% Parameters: Loudspeaker array

%Number of loudspeakers
Nsp = 32;

%Array radius [m]
r_sp = 1.5;

%Loudspeaker interval[rad]
d_sp = 2*pi/Nsp;

%Positions
phi_sp = (0:(Nsp-1))'*d_sp;

x_sp = r_sp*cos(phi_sp);
y_sp = r_sp*sin(phi_sp);
z_sp = zeros(Nsp,1);

%% Parameters: Microphone array

%Number of microphones
Nm = Nsp;

%Array radius [cm]
r_m = 0.25;

%Microphone interval [cm]
d_m = 2*pi/Nm;

%Positions
phi_m = (0:(Nm-1))'*d_m;

x_m = r_m*cos(phi_m);
y_m = r_m*sin(phi_m);
z_m = zeros(Nm,1);

%% Sound sources

%Source type ('ps': point source, 'pw': plane wave)
src_type = 'pw';

%Point source position [m]
rs = 2.5;
theta_s = 90*pi/180;
phi_s = -90*pi/180;

xs = rs*sin(theta_s)*cos(phi_s);
ys = rs*sin(theta_s)*sin(phi_s);
zs = rs*cos(theta_s);

%Plane wave angle [rad]
theta_pw = 90/180*pi;
phi_pw = 90/180*pi;

%Amplitude
amp = 1.0;

%Frequency
freq=1000;

%Wave number
k=2*pi.*freq./c;

%Wave vector of plane wave
kx_pw = k*sin(theta_pw)*cos(phi_pw);
ky_pw = k*sin(theta_pw)*sin(phi_pw);
kz_pw = k*cos(theta_pw);

%% Parameters: Spherical hamonic functions

%Order for WFR filter
Nmax = 16;

size_n = Nmax+1;
size_m = 2*Nmax+1;

n = (0:Nmax).';
m = (-Nmax:Nmax).';

Ynm_mat = zeros(size_n, size_m, Nm);
for iph = 1:Nm
    Ynm_mat(:,:,iph) = sph_harm_mat(Nmax,pi/2,phi_m(iph));
end

%% Received signals of microphone array with rigid spherical buffle

sig = zeros(Nm,1);

%Point source
if strcmp(src_type,'ps')==1
    spec = zeros(size_n,size_m);
    Ynm_ps = sph_harm_mat(Nmax,theta_s,phi_s);
    
    for in=1:size_n
        spec(in,:) = -k./sph_besselh_diff(n(in),1,k*r_m)*sph_besselh(n(in),1,k*rs).*conj(Ynm_ps(in,:));
    end

    for iph=1:Nm
        sig(iph) = sum(sum(spec.*Ynm_mat(1:size_n,1:size_m,iph)));
    end
end

%Plane wave
if strcmp(src_type,'pw')==1
    Ynm_pw = sph_harm_mat(Nmax,theta_pw,phi_pw);
    for iph=1:Nm
        sig(iph) = 0.0;
        for in=1:(Nmax+1)
            sig(iph) = sig(iph) + 4.*pi.*(1i).^(n(in)).*(1i./((k^2).*(r_m^2).*sph_besselh_diff(n(in),1,k*r_m))).*sum(conj(Ynm_pw(in,:)).*Ynm_mat(in,:,iph));
        end
        sig(iph) = amp*sig(iph);
    end
end

%% WFR filtering

%Filter type (1: monopole source assumption, 2: line source assumption)
filter_type = 1;

Fwfr = zeros(Nsp,1);
Fwfr_spec = zeros(size_m,1);

% Reference radius
r_ref = 0.1;

Ynm_0 = Ynm_mat(1:size_n,1:size_m,1);

for im=1:size_m
    coef_nu = 0.0;
    coef_eps = 0.0;
    for in=(abs(m(im))+1):size_n
        coef_nu = coef_nu + (1i^(n(in)+1))./(sph_besselh_diff(n(in),1,k*r_m)).*(Ynm_0(in,n(in)+m(im)+1)^2);
        coef_eps = coef_eps + sph_besselj(n(in),k*r_ref)./(4*pi*r_sp).*sph_besselh(n(in),1,k.*r_sp).*(Ynm_0(in,n(in)+m(im)+1)^2);
    end
    
    if filter_type == 1
        Fwfr_spec(im)=((1i^(m(im)-1))*(r_m^2)).*besselj(m(im),k*r_ref)./(coef_nu.*coef_eps);
    elseif filter_type == 2
        Fwfr_spec(im)=(k*(1i^(m(im)-1)).*(r_m^2))./(pi.*(r_sp/100).*besselh(m(im),k*r_sp))./coef_nu;
        Fwfr_spec(im)=sqrt((2.*pi.*1i)./k).*Fwfr_spec(im);
    end
end

%Windowing circular harmonics
win_m = abs(m)<k*r_m;
Fwfr_spec = Fwfr_spec.*win_m;

for iph=1:Nsp
    Fwfr(iph) = sum(Fwfr_spec.*exp(1i.*m.*phi_sp(iph)));
end

%Convolution
Fwfr_spec = ifft(Fwfr,Nm);
spec_rcv = ifft(sig,Nm);
spec_drv = spec_rcv.*Fwfr_spec;
sig_drv = fft(spec_drv,Nm);

%% Calculate original and reproduced sound field
%Reproduced sound field
dist = zeros(Nx,Ny);

%Original sound field
dist_org = zeros(Nx,Ny);

for iy=1:Ny
    for ix=1:Nx
        %Reproduced sound field
        H_sp2r = green3d(xr(ix), yr(iy), zr, x_sp, y_sp, z_sp, k);
        dist(ix,iy) = sum(H_sp2r.*sig_drv);
        
        %Original sound field
        %Plane wave
        if strcmp(src_type,'ps')==1
            H_s2r = green3d(xr(ix), yr(iy), zr, xs, ys, zs, k);
            dist_org(ix,iy) = sum(amp.*H_s2r);
        elseif strcmp(src_type,'pw')==1
            dist_org(ix,iy) = amp*exp(1i*k*(xr(ix)*sin(theta_pw)*cos(phi_pw)+yr(iy)*sin(theta_pw)*sin(phi_pw)));
        end
    end
end

%Normalize amplitude
dist_org = dist_org./abs(dist_org(Nx/2+1,Ny/2+1));
dist = dist./abs(dist(Nx/2+1,Ny/2+1));

%Normalized error distribution
dist_err = 10*log10(abs(dist - dist_org).^2./abs(dist_org).^2);

%% Draw figures

[XX,YY] = meshgrid(xr, yr);
fig_pos=[0, 0, 600, 500];
fig_prm.label_x = 'x [m]'; fig_prm.label_y = 'y [m]';

%Pressure distribution
fig_prm.zrange=[-4.0,4.0];
fig_dist(figure(1),wloc4_l1(1,:),XX,YY,real(dist_org.'),fig_prm);
fig_dist(figure(2),wloc4_l1(2,:),XX,YY,real(dist.'),fig_prm,x_sp,y_sp);

%Error distribution
fig_prm.zrange=[-30.0,0.0];
fig_contour(figure(3),wloc4_l2(2,:),XX,YY,dist_err.',fig_prm,x_sp,y_sp);

%% Terminate

% rmpath('../lib');
