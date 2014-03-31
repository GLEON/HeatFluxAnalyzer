function [lw,LWo,lwnet] = calc_lwnet(Jday,lat,press,ta,rh,sw,ts)    
    clrSW = clearSkySW(Jday,lat,press,ta,rh*.01);
    clf = 1-sw./clrSW;
    ltI = lt(clf,0);
    clf(ltI) = 0;
    gtI = gt(clf,1);
    clf(gtI) = 1;
    month = datevec(Jday);
    month = month(:,2);
    lw = Mdl_LW(ta,rh*0.01,clf,month);
    Tk = ts + 273.13;     % sample Ts to match LW
    emiss = 0.972;
    S_B = 5.67E-8;
    LWo = S_B*emiss*Tk.^4;
    lwnet = lw-LWo;
end

function clrSW = clearSkySW(time,lat,pressure,temperature,RH)
    if eq(nargin,2)
        pressure = 1013;
        temperature = 20;
        RH       = .7;
    else
        RH = RH*0.01;
    end
    % % resize vectors
    time = reshape(time,length(time),1);
    if gt(length(time),length(lat))
        lat  = ones(length(time),1)*lat(1);
    else
        lat = reshape(lat,length(lat),1);
    end
    if gt(length(time),length(pressure))
        pressure  = ones(length(time),1)*pressure(1);
    else
        pressure = reshape(pressure,length(pressure),1);
    end
    if gt(length(time),length(temperature))
        temperature  = ones(length(time),1)*temperature(1);
    else
        temperature = reshape(temperature,length(temperature),1);
    end
    if gt(length(time),length(RH))
        RH  = ones(length(time),1)*RH(1);
    else
        RH = reshape(RH,length(RH),1);
    end

    % % integrating lat and long into a timezone calculation...
    % URL = ['http://www.askgeo.com/api/334001/jk6n0p5jknsbe5gq0dh1su5bsr/'...
    %     'timezone.json?points=' num2str(lat) '%2C' num2str(long)];
    % str = urlread(URL);
    % currentStI = strfind(str,'"currentOffsetMs"');
    % latitudeStI = strfind(str,',"latitude"');
    % GMT_Off = str2double(str(currentStI+18:latitudeStI)); % in ms, UTC 
    % hr_Off = GMT_Off/3600/1000;     % hour offset from GMT
    % calculate SW from clear sky
    % from Crawford and Duchon 1991 *** key reference here ***
    t_noon = 12.5;          % solar noon (actually will depend on timezone)
    t  = datevec(time);
    n  = floor(time-datenum(t(:,1),1,0));  % Julian day
    t  = t(:,4);            % time (hours)

    cosN = 1+0.034*cos(2*pi*(n-1)/365);
    I0 = 1353*cosN.*cosN;     % Meyers & Dale 1983

    H = (pi/12)*(t_noon-t);   % t_noon is solar noon, t is the local solar time
    % e.g. t = 12.5 and H = 0 at local solar noon

    sigma = 180/pi*(0.006918-.399912*cos(2*pi/365*(n-1))+0.070257*...
        sin(2*pi/365*(n-1))-.006758*cos(4*pi/365*(n-1))+...
        0.000907*sin(4*pi/365*(n-1))-0.002697*cos(6*pi/365*(n-1))+...
        0.00148*sin(6*pi/365*(n-1)));

    % - - cosine of solar zenith angle - -
    sin1 = sin(lat*2*pi/360);
    sin2 = sin(sigma*2*pi/360);
    cos1 = cos(lat*2*pi/360);
    cos2 = cos(sigma*2*pi/360);
    cos3 = cos(H);
    cosZ = sin1.*sin2+cos1.*cos2.*cos3;
    
    p = pressure;                       % in millibars
    m1 = 1224*cosZ.*cosZ+1;
    m = 35*m1.^-.5;                     % optical air mass at p = 1013 mb

    Tr1 = m.*(0.000949*p+0.051);
    T_rT_pg = 1.021-0.084*Tr1.^.5;      % Atwater & Brown 1974

    Es = SatVaporFromTemp(temperature)*1000;
    E  = (RH.*Es);
    Td1 = 243.5*log(E/6.112);
    Td2 = 17.67-log(E/6.112);
    T_d = Td1./Td2;                     % dewpoint (degrees c)  
    T_d = T_d*9/5+32;                   % dewpoint (degrees F)
    %lat=46
    G   = getSmithGamma(lat,time);      % empirical constant dependent upon 
                                        % time of the year and latitude (Smith
                                        % 1966 tale 1
    T_a = m;
    for i = 1:length(m)
        T_a(i) = 0.935^m(i);            % Meyers & Dale 1993
    end               

    mu = exp(0.1133-log(G+1)+0.0393*T_d); % precipitable water

    mu1 = mu.*m;
    T_w = 1-.077*mu1.^0.3;

    ISW = I0.*cosZ.*T_rT_pg.*T_w.*T_a;
    useI = gt(ISW,0);

    clrSW = zeros(length(time),1);
    clrSW(useI) = ISW(useI);
end

function Pressure = SatVaporFromTemp( Temperature )
    %Empirical fit from P.R. Lowe, 1976. An Approximating Polynomial for the
    %Computation of Saturation Vapor Pressure. Journal of Applied Meteorology
    T = Temperature;
    a0 = 6.107799961;
    a1 = 4.436518521E-1;
    a2 = 1.428945805E-2;
    a3 = 2.650648471E-4;
    a4 = 3.031240396E-6;
    a5 = 2.034080948E-8;
    a6 = 6.136820929E-11;
    Pressure = a0 + T.*(a1+T.*(a2+T.*(a3+T.*(a4+T.*(a5*a6*T)))));
    Pressure = Pressure/1000; %now in bar
end

function SmithGamma = getSmithGamma(lat,time)
    if lt(lat,0)
        time = time+182;    % convert for S Hemisphere
        lat = abs(lat);
    end
    yr = datevec(time);
    yr = yr(:,1);
    time = time-datenum(yr,0,1);
    latDeg = 5:10:85;
    season = [-10 81 173 264 355 446];  % winter,spring,summer,fall,winter,spring
    gammaTable = [3.37 2.85 2.8 2.64 3.37 2.85;...
        2.99 3.02 2.7 2.93 2.99 3.02;...
        3.6 3 2.98 2.93 3.6 2.98;...
        3.04 3.11 2.92 2.94 3.04 3.11;...
        2.7 2.95 2.77 2.71 2.7 2.95;...
        2.52 3.07 2.67 2.93 2.52 3.07;...
        1.76 2.69 2.61 2.61 1.76 2.69;...
        1.6 1.67 2.24 2.63 1.6 1.67;...
        1.11 1.44 1.94 2.02 1.11 1.44];
    SmithGamma = interp2(latDeg,season,gammaTable',lat,time);
end
