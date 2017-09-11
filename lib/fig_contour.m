%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Shoichi Koyama (koyama.shoichi@ieee.org) 2016.09.01
%
% Draw pressure distribution
% USAGE "fig_dist(h,X,Y,D,prm,plt_x,plt_y)"
% INPUT
%   h             Figure handler
%   pos           Figure position
%   X, Y          Grid on x and y axes
%   D             Data
%   prm           Parameters 
%   plt_x, plt_y  Plot points (optional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rtn = fig_contour(h,pos,X,Y,D,prm,plt_x,plt_y)

if ~isfield(prm,'level_d'); prm.level_d = 4; end;
if ~isfield(prm,'level_min'); prm.level_min = min(prm.zrange); end;
if ~isfield(prm,'level_max'); prm.level_max = max(prm.zrange); end;
if ~isfield(prm,'plt_linewidth'); prm.plt_linewidth = 2; end;
if ~isfield(prm,'plt_markersize'); prm.plt_markersize = 6; end;
if ~isfield(prm,'plt_markeredgecolor'); prm.plt_markeredgecolor = 'w'; end;
if ~isfield(prm,'plt_markerfacecolor'); prm.plt_markerfacecolor = 'k'; end;

set(h,'Position',pos);
hold on;
contourf(X,Y,D,prm.level_min:prm.level_d:prm.level_max);
if nargin==8
    plot(plt_x,plt_y,'o','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','w','MarkerFaceColor','k');
end
hold off;
caxis(prm.zrange);
axis equal;
axis tight;
colormap(flipud(pink));
colorbar;
xlabel(prm.label_x); ylabel(prm.label_y);

end