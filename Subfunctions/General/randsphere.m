function [vec,x,y,z] = randsphere(sphererad)
% https://www.mathworks.com/help/matlab/math/numbers-placed-randomly-within-volume-of-sphere.html#bt2k246

%Generate random number within sphere
phi = asin(2*rand()-1); %Elevation in open interval of (−π/2,π/2) but not uniformly distributed
theta = 2*pi*rand(); %Azimuth in open interval of (0,π/2), uniformly distributed
rad = sphererad*(rand().^1/3); %Sphere radii not uniformly distributed

%Convert to x, y, z coordinates
[x,y,z] = sph2cart(theta,phi,rad);
vec = [x;y;z];
end