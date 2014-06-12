function LW = Mdl_LW(airT,RH,clf,month)
    S_B = 5.67E-8;
    vp = VaporPressure(airT,RH);   % air vapor pressure (mbars)
    T_k = 273.13+airT;    
    cl1 = 1-clf;
    cl2 = 1.22+0.06*sin((month+2)*pi/6);
    cl3 = vp./T_k;
    T_k1 = T_k.^4;
    LW = T_k1.*(clf+cl1.*cl2.*cl3.^(1/7))*S_B;
end

function [VaporPressure ] = VaporPressure( Temperature,RelativeHumidity )
    %Temperature in °C, relative Humidity as decimal (1 > Rh > 0)
    %VaporPressure in mb
    satPressure = SatVaporFromTemp( Temperature )*1000; %satP in mb
    VaporPressure = RelativeHumidity.*satPressure;
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
    Pressure = a0 + T.*(a1+T.*(a2+T.*(a3+T.*(a4+T.*(a5+a6*T)))));
    Pressure = Pressure/1000; %now in bar
end
