function sw_alb = sw_albedo(Jday,lat)
% INPUTS:
%   lat: latitude in degrees North.
%   Jday: Julian Date; serial date number which represents the whole and 
%         fractional number of days from January 0, 0000.
%
% OUTPUTS:
%   sw_alb: albedo for short wave radiation for given latitude and time of year

dateV = datevec(Jday); % date vector
doy = Jday - datenum(dateV(:,1),0,0);

% Assign Parameter Values %
Tropic = 23.45; 
YrDays = 365.0;       
DclDay = 284.0;
DgCrcl = 360.0; 
DgToRd = 57.29577951; 
RefInd = 1.34; % This is a mild function of T and S. http://scubageek.com/articles/wwwh2o.html

Lat = lat/DgToRd;                                  % Start year at 0.0
z = DclDay + floor(doy);
x = DgCrcl*z/YrDays/DgToRd;
y = sin(x);
Decl = Tropic/DgToRd*y;

% Hour-angle calculation where time is hours from midnight  %
HrAng = ((doy-floor(doy))*24 -12)*15.0/DgToRd;

% Zenith angle calculation %
Zenith = acos(sin(Decl)*sin(Lat)+cos(Decl)*cos(Lat).*cos(HrAng));

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