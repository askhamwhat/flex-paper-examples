
function [r,d,d2] = starfish3arms(t)
%STARFISH return position, first and second derivatives of a starfish 
% domain with the parameterization 
%
% x(t) = (1+amp*cos(narms*(t+phi)))*cos(t)
% y(t) = (1+amp*cos(narms*(t+phi)))*sin(t)
%
% Syntax: [r,d,d2] = starfish(t,narms,amp,ctr,phi,scale)
%
% Input:
%   t - array of points (in [0,2pi])
%
% Optional input:
%   narms - integer, number of arms on starfish (5)
%   amp - float, amplitude of starfish arms relative to radius of length 1
%               (0.3)
%   ctr - float(2), x0,y0 coordinates of center of starfish ( [0,0] )
%   phi - float, phase shift (0)
%   scale - scaling factor (1.0)
%
% Output:
%   r - 2 x numel(t) array of positions, r(:,i) = [x(t(i)); y(t(i))]
%   d - 2 x numel(t) array of t derivative of r 
%   d2 - 2 x numel(t) array of second t derivative of r 
%
% Examples:
%   [r,d,d2] = starfish(t); % get default settings 
%   [r,d,d2] = starfish(t,narms,[],ctr,[],scale); % change some settings


narms = 3;
amp = 0.3;

ct = cos(t);
st = sin(t);
cnt = cos(narms*(t));
snt = sin(narms*(t));

xs = (1+amp*cnt).*ct;
ys = (1+amp*cnt).*st;
dxs = -(1+amp*cnt).*st-narms*amp*snt.*ct;
dys = (1+amp*cnt).*ct-narms*amp*snt.*st;

d2xs = -dys-narms*amp*(narms*cnt.*ct-snt.*st);
d2ys = dxs-narms*amp*(narms*cnt.*st+snt.*ct);

r = [(xs(:)).'; (ys(:)).'];
d = [(dxs(:)).'; (dys(:)).'];
d2 = [(d2xs(:)).'; (d2ys(:)).'];

end

