function [kx,ky,NAt]= k_vector(xi,yi,H,LEDp,nglass,t,theta,xint,yint,total)
% Calculate k vector and illumination NA corresponding to each LED
% Input:
%       xi,yi: relative coordinate of lit LEDs
%       xint,yint: offset of initial LED to the patch center, in mm
%       (x is along horizonal direction; y is along vertical direction)
%       H: distance between LEDs and sample, in mm
%       LEDp: distance between adjacent LEDs, in mm
%       nglass: refraction index of glass substrate
%       t: glass thickness, in mm
%       theta: rotation angle of LED array to the camera sensor frame, in degree
%       total: total number of lit LEDs
% Output:
%       kx,ky: normalized wavevectors, which should times k0 later
%       NAt: illumination NA of each LED

kx  = zeros(1,total);
ky  = kx;
NAt = kx;

for tt = 1:total
    x0 = xint+xi(1,tt)*LEDp; % from rightmost postion
    y0 = yint+yi(1,tt)*LEDp; % from topmost postion
    x1 = x0*cos(theta*pi/180)-y0*sin(theta*pi/180);
    y1 = x0*sin(theta*pi/180)+y0*cos(theta*pi/180);
    [kx(1,tt),ky(1,tt),NAt(1,tt)] = calculate(x1,y1,H,t,nglass);
end

end

function [kx,ky,NAt]=calculate(x0,y0,H,h,n)
% iterative root finder for LED offset to account for the effect of microslide 
% this function converges in one iteration for H=83, h=1, n=1.45 & x0,y0 up to 83.

l = sqrt(x0^2+y0^2); % distance of LED from origin
thetal = atan2(y0,x0); % angle of LED in x-y plane

xoff = 0; % initial guess where beam enters bottom of slide
thetag = -asin(l/sqrt(l^2+H^2)/n); % get angle of beam in glass from Snell's law
xint = h*tan(thetag); % find where the beam exits the top of the slide
xoff = xoff-xint; % modify guess where beam enters bottom of slide by this amount

% repeat the above procedure until the beam exits the top of the slide
% within 1 micron of center
while abs(xint) > .001
    thetag = -asin((l-xoff)/sqrt((l-xoff)^2+H^2)/n);
    xint = xoff+h*tan(thetag);
    xoff = xoff-xint;
end

% angle under the glass and angle over the cover slip
% FPM treats this as the angle in the sample so assumes the sample has
% refractive index = 1.0
theta = asin((l-xoff)/sqrt((l-xoff)^2+H^2));

NAt = abs(sin(theta));
kx  = -NAt*cos(thetal);
ky  = -NAt*sin(thetal);
end