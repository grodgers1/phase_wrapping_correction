function [m_dpc_sino] = make_dpc_sino(m_sino,delta,d_t,p2)
%MAKE_DPC_SINO takes a sinogram of projected thickness and turns it into
% sinogram with interference pattern phase shift

%   m_sino: sinogram of modelled projected thickness (from 'model_cyl_sino.m')
%   delta: modelled delta value
%   d_t: talbot distance [m]
%   p2: grating period [m]

%   m_dpc_sino: modeled dpc sino, giving interference pattern phase shift

m_rad_sino = m_sino * 2 * pi * delta; 
temp = interp1(1:size(m_rad_sino,1), m_rad_sino, 1:0.5:size(m_rad_sino,1)); % expand vector
temp2 = circshift(temp,2)-temp; % take derivative
m_diff_sino = interp1(1:size(temp2,1), temp2, 1:2:size(temp2,1)); % reduce to original size
m_dpc_sino = -(d_t / p2) * m_diff_sino; % modeled dpc signal

end

