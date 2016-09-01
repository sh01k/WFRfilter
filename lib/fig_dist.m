%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Shoichi Koyama (koyama.shoichi@ieee.org) 2016.09.01
%
% Draw pressure distribution
% USAGE "fig_dist(h,pos,X,Y,D,prm,plt_x,plt_y)"
% INPUT
%   h             Figure handler
%   pos           Figure position
%   X, Y          Grid on x and y axes
%   D             Data
%   prm           Parameters 
%   plt_x, plt_y  Plot points (optional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rtn = fig_dist(h,pos,X,Y,D,prm,plt_x,plt_y)

if ~isfield(prm,'plt_linewidth'); prm.plt_linewidth = 2; end;
if ~isfield(prm,'plt_markersize'); prm.plt_markersize = 6; end;
if ~isfield(prm,'plt_markeredgecolor'); prm.plt_markeredgecolor = 'w'; end;
if ~isfield(prm,'plt_markerfacecolor'); prm.plt_markerfacecolor = 'k'; end;

set(h,'Position',pos,'PaperUnits','centimeters','PaperPosition',[1 1 10 10]);
hold on;
surf(X,Y,D);
if nargin==8
    plot3(plt_x,plt_y,ones(length(plt_x),1)*max(max(abs(D))),'o','LineWidth',prm.plt_linewidth,'MarkerSize',prm.plt_markersize,'MarkerEdgeColor',prm.plt_markeredgecolor,'MarkerFaceColor',prm.plt_markerfacecolor);
end
hold off;
shading flat;
view(0,90);
caxis(prm.zrange);
axis equal;
axis tight;
colormap(flipud(pink));
colorbar;
xlabel(prm.label_x); ylabel(prm.label_y);

end