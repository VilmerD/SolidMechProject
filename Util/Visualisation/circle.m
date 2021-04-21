function circle(p, r, linestyle)
% Plots a circle in 2d-space
%
% Inputs
%   p           mid point of circle [xc, yc]
%
%   r           radius of circle
%
%   linestyle   linestyle of circle, optional
if nargin < 3
   linestyle = 'k'; 
end
tt = linspace(0, 2*pi);
plot(p(1) + r*cos(tt), p(2) + r*sin(tt), linestyle);
end