function [l_edge,r_edge,l_pw,r_pw,r,x0,gof_vec] = find_edges(edge_sino,pw_sino,bump,smoothedges,smoothfactor,usefits)
%FIND_EDGES Finds the edges of the sinogram and the extent of phase
%wrapping. Bump is an extra buffer from the edge of the extent of phase
%wrapping.

%   edge_sino: sino to be used to find the geometry of the specimen (later these are
%used in the modeling). 
% This is usually the abs sino, but it could be the phase sino plus some
% bump. Or it could be the derivative of the abs sino. The current script
% finds the max edge_sino on the left and right sides
%   pw_sino: the phase wrapped sino to find extent of wrapping
%   bump : pixel bump when finding extent of phase wrapping (often max of
%pw_sino is a few pixels away from where the wrapping begins)
%   smoothedges: 1=on, 0=off
%   smoothfactor: Over what window smoothing ('rlowess') is done.
%   usefits: 1=on, 0=off. If on, this does a fit of the found edges and
%uses the fit values for modeling the container.

% % This version was made on May 31, 2018 by cleaning up the old script

l_pw = zeros(size(edge_sino,2),1); % first phase wrapped pixel on left
r_pw = zeros(size(edge_sino,2),1); % first phase wrapped pixel on right
l_edge = zeros(size(edge_sino,2),1); % left edge of sample position
r_edge = zeros(size(edge_sino,2),1); % right edge of sample position

% find the max of the temp_sino on the left and right side of each line
for ang = 1:size(edge_sino, 2) % loop over angles
    l_edge(ang) = find(edge_sino(1:round(end/2),ang) == ...
        max(edge_sino(1:round(end/2),ang)),1);
    r_edge(ang) = find(edge_sino(round(end/2):end,ang) == ...
        max(edge_sino(round(end/2):end,ang)),1,'last')+round(size(edge_sino,1)/2)-1;
    l_pw(ang) = find(pw_sino(1:round(end/2),ang) == ...
        max(pw_sino(1:round(end/2),ang)),1) + bump;
    r_pw(ang) = find(pw_sino(round(end/2):end,ang) == ...
        max(pw_sino(round(end/2):end,ang)),1,'last')+round(size(pw_sino,1)/2)-1 - bump;
end

% optionally smooth these edges
if smoothedges == 1
    l_edge = round(smooth(l_edge, smoothfactor, 'rlowess'));
    r_edge = round(smooth(r_edge, smoothfactor, 'rlowess'));
    l_pw = round(smooth(l_pw, smoothfactor, 'rlowess'));
    r_pw = round(smooth(r_pw, smoothfactor, 'rlowess'));
end

% find center and radius
x0 = size(edge_sino, 1)/2 - (r_edge + l_edge)/2; % relative to middle of projection
r = mean((r_edge - l_edge)/2);

% Fit the edges
xdata = (1:size(edge_sino, 2))';
[f_ledge, gof_ledge] = fit(xdata, l_edge, 'fourier1','Robust','on');
[f_redge, gof_redge] = fit(xdata, r_edge, 'fourier1','Robust','on');
[f_lpw, gof_lpw] = fit(xdata, l_pw, 'fourier1','Robust','on');
[f_rpw, gof_rpw] = fit(xdata, r_pw, 'fourier1','Robust','on');
gof_vec = [gof_ledge.rsquare, gof_redge.rsquare, gof_lpw.rsquare, gof_rpw.rsquare];

if usefits == 1
    l_edge = round(feval(f_ledge,xdata));
    r_edge = round(feval(f_redge,xdata));
    l_pw = round(feval(f_lpw,xdata));
    r_pw = round(feval(f_rpw,xdata));
    x0 = size(edge_sino, 1)/2 - (r_edge + l_edge)/2; % relative to middle of projection
    r = mean((r_edge - l_edge)/2);
end

end

