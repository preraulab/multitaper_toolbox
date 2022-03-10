function c = colorbar_noresize(varargin)
%COLORBAR_NORESIZE Makes a colorbar that does not resize the axis
%   colorbar_noresize(ax,<colorbar arguments>)
%
% Copyright Michael Prerau 2018
if nargin==0 || ~strcmpi(class(varargin{1}),'matlab.graphics.axis.Axes')
    ax=gca;
else
    ax=varargin{1};
end

pos = ax.Position;
c = colorbar(ax, varargin{2:end});
ax.Position = pos;

end

