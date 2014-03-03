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

if lat == 0
    sw_alb = 0.08;
end

if lat > 0;
    sw_alb = 0.08 + 0.02.*sin(((2.*pi./365).*doy) + (pi./2));
end

if lat < 0
    sw_alb = 0.08 - 0.02.*sin(((2.*pi./365).*doy) + (pi./2));
end
            
end