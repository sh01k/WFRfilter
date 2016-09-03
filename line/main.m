%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Shoichi Koyama (koyama.shoichi@ieee.org) 2016.09.01
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

%Sound velocity [m/s]
c=340.29;

%% Parameters: Reproduced area

%Size of reproduced area [m]
lenX=4.2;
lenY=4.2;

%Interval of reproduced area [m]
dx=0.015;
dy=0.015;

%Number of samples
Nx=round(lenX/dx);
Ny=round(lenY/dy);

%Positions
xr = (((0:Nx-1)-Nx/2).*dx)';
yr = ((0:Ny-1).*dy-0.2)';
zr = 0;

%% Parameters: Loudspeaker array

%Number of loudspeakers
Nsp = 64;

%Loudspeaker interval [m]
d_sp = 0.06;

%Array size [m]
len_sp = Nsp*d_sp;

%Positions
if mod(Nsp,2)==0
    x_sp = (((0:Nsp-1)-Nsp/2+0.005).*d_sp)';
else
    x_sp = (((0:Nsp-1)-(Nsp-1)/2).*d_sp)';
end
y_sp = zeros(Nsp,1);
z_sp = zeros(Nsp,1);

%% Parameters: Microphone array

%Number of microphones
Nm = Nsp;

%Microphone interval [m]
d_m = d_sp;

%Array size [m]
len_m = len_sp;

%Positions
x_m = x_sp;
y_m = y_sp;
z_m = z_sp;

%% Sound sources

%Number of sources
Nsrc = 1;

%Source type ('ps': point source, 'pw': plane wave)
src_type = 'ps';

%Source position [m]
xs = -0.4;
ys = -1.0;
zs = 0.0;

%Plane wave angle [rad]
theta_pw = 5/180*pi;

%Amplitude
amp = 10.0;

%Frequency
freq=1000;

%Wave number
k=2*pi*freq/c;

%% Spatial window function

win_sp=ones(Nsp,1);

taper_point=ceil(Nsp*0.10);
taper_vec=0.5-0.5.*cos(pi.*(1:taper_point)./taper_point);

win_sp(1:taper_point)=taper_vec;
win_sp(Nsp-taper_point+1:Nsp)=taper_vec(taper_point:-1:1);

%% Received signals of microphone array

sig = zeros(Nm,1);

%Point source
if strcmp(src_type,'ps')==1
    for iis=1:Nsrc
        H_s2m = green3d(x_m,y_m,z_m,xs(iis),ys(iis),zs(iis),k);
        sig = sig + amp(iis).*H_s2m;
    end
end
%Plane wave
if strcmp(src_type,'pw')==1
    for iis=1:Nsrc
        sig = sig + amp(iis).*exp(1i*(k*sin(theta_pw(iis))*x_m+k*cos(theta_pw(iis))*y_m));
    end
end

%% WFR filtering

%Spatial filter length
wfr_filter_len = 128;

%Spatial DFT length
Nf = 256;

%Reference line [m]
y_ref = 2.0;

%Shift [m]
shift.x = 0.0;
shift.y = 0.0;

%Spatial frequency
if mod(Nf,2)==0
    kx = (2*pi*([0:Nf/2-1, -Nf/2:-1]/Nf)/d_m)';
else
    kx = (2*pi*([0:(Nf-1)/2, -(Nf-1)/2:-1]/Nf)/d_m)';
end

%Spatial frequency window function
kx_c = k;
win_kx = (abs(kx)<kx_c);

%WFR filter
ky=sqrt(k^2-kx.^2);
Fwfr_spec=-4i.*exp(1i.*ky.*y_ref)./besselh(0,1,ky.*y_ref).*exp(-1i.*(kx.*shift.x+ky.*shift.y));

%Windowing in spatio-temporal frequency domain
Fwfr_spec = Fwfr_spec.*win_kx;

%Padding zero
Fwfr = fft(Fwfr_spec,Nf);
Fwfr(1+wfr_filter_len/2:Nf-wfr_filter_len/2) = zeros(1,Nf-wfr_filter_len);

%Filter convolution
Fwfr_spec = ifft(Fwfr, Nf);

sig_rcv = zeros(Nf,1);
sig_rcv((Nf-Nm)/2+1:Nf-(Nf-Nm)/2)=sig;
sig_rcv = [sig_rcv(Nf/2+1:Nf); sig_rcv(1:Nf/2)];

spec_rcv = ifft(sig_rcv, Nf);
spec_drv = spec_rcv.*Fwfr_spec;
sig_drv = fft(spec_drv, Nf);

%Permutate
sig_drv = [sig_drv(Nf/2+1:Nf); sig_drv(1:Nf/2)];
sig_drv = sig_drv((Nf-Nm)/2+1:Nf-(Nf-Nm)/2);

%Tapering window
sig_drv = sig_drv.*win_sp;

%% Calculate original and reproduced sound field

%Reproduced sound field
dist = zeros(Nx,Ny);

%Original sound field
dist_org = zeros(Nx,Ny);

for iy=1:Ny
    for ix=1:Nx
        %Reproduced sound field
        H_sp2r = green3d(xr(ix),yr(iy),zr,x_sp,y_sp,z_sp,k);
        dist(ix,iy) = sum(H_sp2r.*sig_drv);
        
        %Original sound field
        if strcmp(src_type,'ps')==1
            H_s2r = green3d(xr(ix),yr(iy),zr,xs+shift.x,ys+shift.y,zs,k);
            dist_org(ix,iy) = sum(amp.*H_s2r);
            
        elseif strcmp(src_type,'pw')==1
            dist_org(ix,iy) = sum(amp.*exp(1i*(k*sin(theta_pw)*xr(ix)+k*cos(theta_pw)*yr(iy))));
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
