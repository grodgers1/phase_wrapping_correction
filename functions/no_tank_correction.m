function [unwrap_sino, gof_vec,r,x0] = no_tank_correction(wrap_sino, abs_sino, d_t, p2, delta, buffer, bump, smoothedges, smoothfactor,usefits)
%NO_TANK_CORRECTION Corrects for phase wrapping at edge of sino by modeling
% a cylindrical sample and replacing wrapped regions
%   Inputs:
%   wrap_sino: the phase wrapped sino to correct
%   abs_sino: the accompanying absorption sinogram for finding edges
%   d_t: Talbot distance [m]
%   p2: period of 2nd grating (or interference fringes from g1) [m]
%   delta: decrement of the real part of the refractive index for modelling
%   buffer: around the edges and wrapping extent for definition of
%replacement window.
%   bump: in case adjustments needed for finding phase wrapping. (see
%'find_edges.m')

%   Outputs:
%   unwrap_sino: corrected sinogram
%   gof_vec: gives a goodness of fit vector, describing the quality of the
%fit done in 'find_edges.m' for l_edge,r_edge,l_pw,r_pw. Bigger gof_vec
%means the fit was better, worse usually means it was noisy/not cylindrical
%   r: radius of cylinder
%   x0: center position of cylinder for each projection angle

[npix, nproj] = size(wrap_sino);
%% Finding sample center, edges, and radius in sinogram 
[l_edge,r_edge,l_pw,r_pw,r,x0,gof_vec] = find_edges(abs_sino,wrap_sino,bump,smoothedges,smoothfactor,usefits);
%% Model of a uniform, cylindrical sample
m_sino = model_cyl_sino(x0,r,npix,nproj);
%% Make modeled DPC sinogram
[m_dpc_sino] = make_dpc_sino(m_sino,delta,d_t,p2);
%% Replace phase wrapped data with modeled data
unwrap_sino = wrap_sino; % initialize corrected sinogram
for p = 1:size(wrap_sino,2)
unwrap_sino((l_edge(p)-buffer):(l_pw(p)+buffer),p) = m_dpc_sino((l_edge(p)-buffer):(l_pw(p)+buffer),p); 
unwrap_sino((r_pw(p)-buffer):(r_edge(p)+buffer),p) = m_dpc_sino((r_pw(p)-buffer):(r_edge(p)+buffer),p);
end

end

