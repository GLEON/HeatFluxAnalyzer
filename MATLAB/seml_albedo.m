function [Albdo] = seml_albedo(doy,Lt)
% This program calculates the albedo on a flat lake from Fresnal's
% Law (Neumann & Pierson, 1966).  The angle of refraction is given
% by Snell's Law.  The index of refraction is assumed to be
% independent of wavelength.  The Zenith angle is calculated from
% Milankovich (1930).

%This is copied straight from seml by FRAM 1/6/06

% Assign Parameter Values %
Tropic = 23.45; 
YrDays = 365.0;       
DclDay = 284.0;
DgCrcl = 360.0; 
DgToRd = 57.29577951; 
RefInd = 1.34; % This is a mild function of T and S. http://scubageek.com/articles/wwwh2o.html

Lat = Lt/DgToRd;                                  % Start year at 0.0
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
Albdo = 0.5 * (A1./A2 + A3./A4);

% Set albedo to 1 if greater than 1 %
Albdo(Albdo>1) = 1;