%% name: wrap
%%
%% syntax: [datout] = wrap(datin,wdim)
%%
%% description:
%%		imgin/out - input/output 2D data set
%%		wdim - wrap along which dimension

function [datout] = wrap(datin,wdim);

if (nargin ~= 2)
    fprintf('Usage:\n');
    fprintf('[datout] = wrap(datin,wrapdimension);\n');
    return;
end

datout = ((datin+pi) - 2*pi*floor((datin+pi)/(2*pi))) - pi;
