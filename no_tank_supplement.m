%% Supplementary Material: Program Code



%% Main function for phase wrapping correction:

function [unwrap_sino,gof_vec,r,x0] = phase_wrapping_correction(wrap_sino, abs_sino, d_t, p2, delta, buffer, bump, smoothedges, smoothfactor, usefits)
%   PHASE_WRAPPING_CORRECTION Corrects phase wrapping at edge of sino by
%   modeling a cylindrical sample and replacing wrapped regions
%   Inputs: wrap_sino: the phase wrapped sino to correct 
%   abs_sino: accompanying absorption sinogram for finding edges 
%   d_t: Talbot distance [m]
%   p2: period of 2nd grating (or interference fringes from g1) [m]
%   delta: decrement of the real part of the refractive index 
%   for modelling 
%   buffer: around the edges and wrapping extent for definition 
%       of replacement window.
%   bump: in case adjustments needed for finding phase wrapping. (see
%   'find_edges.m')

%   Outputs: unwrap_sino: corrected sinogram gof_vec: gives a 
%   goodness of fit vector, describing the quality of the fit
%   done in 'find_edges.m' for l_edge,r_edge,l_pw,r_pw.
%   Bigger gof_vec means the fit was better, worse usually 
%   means it was noisy/not cylindrical
%   r: radius of cylinder 
%   x0: center position of cylinder for each projection angle 

    [npix,nproj] = size(wrap_sino);
    % Finding sample center, edges, and radius in sinogram
    [l_edge,r_edge,l_pw,r_pw,r,x0,gof_vec] = find_edges(abs_sino, wrap_sino, bump, smoothedges, smoothfactor, usefits);
    % Model of a uniform, cylindrical sample
    m_sino = model_cyl_sino(x0,r,npix,nproj);
    % Make modeled DPC sinogram
    [m_dpc_sino] = make_dpc_sino(m_sino,delta,d_t,p2);
    % Replace phase wrapped data with modeled data
    unwrap_sino = wrap_sino; % initialize corrected sinogram
    for p = 1:size(wrap_sino,2)
        unwrap_sino((l_edge(p)-buffer):(l_pw(p)+buffer),p) = m_dpc_sino((l_edge(p)-buffer):(l_pw(p)+buffer),p);
        unwrap_sino((r_pw(p)-buffer):(r_edge(p)+buffer),p) = m_dpc_sino((r_pw(p)-buffer):(r_edge(p)+buffer),p);
    end

end


%% Sub-functions called by "phase_wrapping_correction.m"

function [l_edge,r_edge,l_pw,r_pw,r,x0,gof_vec] = find_edges(edge_sino, pw_sino, bump, smoothedges, smoothfactor, usefits)
%FIND_EDGES Finds the edges of the sinogram and the extent of 
%   phase wrapping. Bump is an extra buffer from the edge of
%   the extent of phase wrapping.

%   edge_sino: sino to be used to find the geometry of the specimen 
%   (later these are used in the modeling). 
% This is usually the abs sino, but it could be the phase sino plus 
% some bump. Or it could be the derivative of the abs sino. The 
% current script finds the max edge_sino on the left and right sides
%   pw_sino: the phase wrapped sino to find extent of wrapping
%   bump : pixel bump when finding extent of phase wrapping (often max 
%   of pw_sino is a few pixels away from where the wrapping begins)
%   smoothedges: 1=on, 0=off
%   smoothfactor: Over what window smoothing ('rlowess') is done.
%   usefits: 1=on, 0=off. If on, it does a fit of the found edges and
%   uses the fit values for modeling the container.

    l_pw = zeros(size(sino1,2),1); % first phase wrapped pixel on left
    r_pw = zeros(size(sino1,2),1); % first phase wrapped pixel on right
    l_edge = zeros(size(sino1,2),1); % left edge of sample position
    r_edge = zeros(size(sino1,2),1); % right edge of sample position

    % find the max of the temp_sino on the left and right side of each line
    for ang = 1:size(edge_sino, 2) % loop over angles
        l_edge(ang) = find(edge_sino(1:round(end/2),ang) == max(edge_sino(1:round(end/2),ang)),1);
        r_edge(ang) = find(edge_sino(round(end/2):end,ang) == max(edge_sino(round(end/2):end,ang)),1,'last')+round(size(edge_sino,1)/2)-1;
        l_pw(ang) = find(pw_sino(1:round(end/2),ang) == max(pw_sino(1:round(end/2),ang)),1) + bump;
        r_pw(ang) = find(pw_sino(round(end/2):end,ang) == max(pw_sino(round(end/2):end,ang)),1,'last')+round(size(pw_sino,1)/2)-1 - bump;
    end

    % optionally smooth these edges
    if smoothedges == 1
        l_edge = round(smooth(l_edge, smoothfactor, 'rlowess'));
        r_edge = round(smooth(r_edge, smoothfactor, 'rlowess'));
        l_pw = round(smooth(l_pw, smoothfactor, 'rlowess'));
        r_pw = round(smooth(r_pw, smoothfactor, 'rlowess'));
    end

    % find center and radius
    x0 = size(sino1, 1)/2 - (r_edge + l_edge)/2; % relative to middle of projection
    r = mean((r_edge - l_edge)/2);

    % Fit the edges
    xdata = (1:size(sino1, 2))';
    [f_ledge, gof_ledge] = fit(xdata, l_edge, 'fourier1','Robust','on');
    [f_redge, gof_redge] = fit(xdata, r_edge, 'fourier1','Robust','on');
    [f_lpw, gof_lpw] = fit(xdata, l_pw, 'fourier1','Robust','on');
    [f_rpw, gof_rpw] = fit(xdata, r_pw, 'fourier1','Robust','on');
    gof_vec = [gof_ledge.rsquare, gof_redge.rsquare, gof_lpw.rsquare,gof_rpw.rsquare];

    if usefits == 1
        l_edge = round(feval(f_ledge,xdata));
        r_edge = round(feval(f_redge,xdata));
        l_pw = round(feval(f_lpw,xdata));
        r_pw = round(feval(f_rpw,xdata));
        x0 = size(sino1, 1)/2 - (r_edge + l_edge)/2; % relative to middle of projection
        r = mean((r_edge - l_edge)/2);
    end

end


function [m_sino] = model_cyl_sino(x0, r, npix, nproj)
%MODEL_CYL_SINO Generates a sinogram of proj. thickness of a circle
%with specified radius and center
%   x0: center position at each angle (length(x0) should be nproj)
%   r: radius of cylinder
%   npix: number of detector pixels
%   nproj: number of projections
%   m_sino: will be a sinogram with projected thickness of a cylinder

    uvec = 0:(npix-1);
    uvec = uvec - round(size(uvec,2)/2) + 0.5;
    t = zeros(size(uvec)); % will contain projected thickness of circle 
    m_sino = zeros(npix,nproj);
    for j = 1:nproj
        for i = 1:size(uvec,2)
            if abs(uvec(i)+x0(j)) > r % rays missing model circle
                t(i) = 0;
            else % rays going through model circle
                t(i) = 2 * sqrt(r^2 -(uvec(i)+x0(j))^2);
            end
        end
        m_sino(:,j) = t;
    end

end

function [m_dpc_sino] = make_dpc_sino(m_sino, delta, d_t, p2)
%MAKE_DPC_SINO takes a sinogram of projected thickness and makes a
% sinogram with interference pattern phase shift

%   m_sino: sinogram of modelled projected thickness 
%       (from 'model_cyl_sino.m')
%   delta: modelled delta value
%   d_t: talbot distance [m]
%   p2: grating period [m]
%   m_dpc_sino: modeled dpc sino, giving int. pattern phase shift

    m_rad_sino = m_sino * 2 * pi * delta; 
    temp = interp1(1:size(m_rad_sino,1),m_rad_sino,1:0.5:size(m_rad_sino,1)); % expand vector
    temp2 = circshift(temp,2)-temp; % take derivative
    m_diff_sino = interp1(1:size(temp2,1), temp2, 1:2:size(temp2,1)); 
    % reduce to original size
    m_dpc_sino = -(d_t / p2) * m_diff_sino; % modeled dpc signal
end







