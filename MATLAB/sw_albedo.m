function sw_alb = sw_albedo(Jday,lat)
% INPUTS:
%   lat: latitude in degrees North.
%   Jday: Julian Date; serial date number which represents the whole and 
%         fractional number of days from January 0, 0000.
%
% OUTPUTS:
%   sw_alb: albedo for short wave radiation for given latitude and time of year

Jday = Jday(:);
dateV = datevec(Jday); % date vector
doy = Jday - datenum(dateV(:,1),0,0);
if length(lat) ~= 1;
    lat = lat(1);
end

% Assign Parameter Values %
RefInd = 1.33;               
Decl = 180/pi*(0.006918-.399912*cos(2*pi/365*(doy-1))+0.070257*...
    sin(2*pi/365*(doy-1))-.006758*cos(4*pi/365*(doy-1))+...
    0.000907*sin(4*pi/365*(doy-1))-0.002697*cos(6*pi/365*(doy-1))+...
    0.00148*sin(6*pi/365*(doy-1)));

% Hour-angle calculation where time is hours from midnight  %
t_noon = 12.5; 
t = dateV(:,4);
HrAng = (pi/12)*(t_noon-t); 

% Zenith angle calculation % lat is converted to radians here
sin1 = sin(lat*2*pi/360);
sin2 = sin(Decl*2*pi/360);
cos1 = cos(lat*2*pi/360);
cos2 = cos(Decl*2*pi/360);
cos3 = cos(HrAng);
cosZ = sin1.*sin2+cos1.*cos2.*cos3;
Zenith = acos(cosZ);

% Angle of Refraction calculation based on Snell's Law %
RefAng = asin(sin(Zenith)/RefInd); % Angle of refraction

% Albedo Calculation %
A1 = tan(Zenith-RefAng).^2;
A2 = tan(Zenith+RefAng).^2;
A3 = sin(Zenith-RefAng).^2;
A4 = sin(Zenith+RefAng).^2;
sw_alb = 0.5 * (A1./A2 + A3./A4);

% Set albedo to 1 if greater than 1 %
sw_alb(sw_alb > 1) = 1;

end
            