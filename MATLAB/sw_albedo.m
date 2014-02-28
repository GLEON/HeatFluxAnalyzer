function sw_alb = sw_albedo(Jday,lat)
% INPUTS:
%   lat: latitude in degrees North.
%   Jday: Julian Date; serial date number which represents the whole and 
%         fractional number of days from January 0, 0000.
%
% OUTPUTS:
%   sw_alb: albedo for short wave radiation for given latitude and time of year
%
% Cogley, J.G. The albedo of water as a function of latitude. 1979.
% Monthly Weather Review, 107: 775-781.

dateV = datevec(Jday); % date vector
Month = dateV(:,2);
doy = Jday - datenum(dateV(:,1),0,0);

lat = abs(lat);

% calculate albedo following Cogley (1979)          
albedo = [90,nan,30.1,33.3,25.3,16.7,13.3,15,22.6,31.7,30.1,NaN,NaN,18.1;...
    80,30.1,33.7,26.6,17.8,13.8,12.3,13.2,16.3,23.8,32.9,30.1,NaN,15.7;...
    70,34,28.1,18.5,12.2,9.9,9.5,9.7,11.3,16.3,25.4,33.6,32.5,12.8;...
    60,26.3,19.3,12.7,9.3,8,7.7,7.9,8.8,11.4,17.4,24.9,29.4,11.1;...
    50,17.8,13.1,9.5,7.7,7.1,7,7,7.5,8.8,12,16.9,19.8,9.4;...
    40,12.1,9.7,7.8,6.9,6.6,6.5,6.6,6.8,7.5,9.1,11.6,13.2,8;...
    30,9.1,7.9,7,6.5,6.4,6.4,6.4,6.4,6.8,7.6,8.9,9.7,7.2;...
    20,7.6,7,6.5,6.3,6.3,6.4,6.3,6.3,6.4,6.9,7.5,7.9,6.7;...
    10,6.9,6.5,6.3,6.3,6.5,6.6,6.5,6.3,6.3,6.5,6.8,7,6.5]; % table 6
            
% find the closest latitude band
closest = abs(albedo(:,1) - lat);
[~,Indx] = min(closest');
            
% define the albedo for the latitude defined
Rng = albedo(Indx,2:end-1);
            
% return time series of albedo values 
sw_alb = Rng(Month);
sw_alb = sw_alb./100; % convert to fraction 
sw_alb = sw_alb(:); % ensure column vector
            
% if there are any missing values, calculate with other method
idx = isnan(sw_alb);

if sum(idx) > 0;
    
    % northern hemisphere
    if lat > 0;
        sw_alb(idx) = 0.08 + 0.02.*sin((((2.*pi)./365).*doy(idx)) + (pi./2));
    end
    % equator
    if lat == 0;
        sw_alb(isnan(sw_alb)) = 0.08;
    end
    % southern hemisphere
    if lat < 0;
        sw_alb(idx) = 0.08 - 0.02.*sin((((2.*pi)./365).*doy(idx)) + (pi./2));
    end

end
end