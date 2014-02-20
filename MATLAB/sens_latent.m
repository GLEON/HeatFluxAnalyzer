function mm = sens_latent(ts,Uz,ta,rh,hu,ht,hq,alt)

% ts: surface temperature, degrees C.
% Uz: wind speed, m/s.
% ta: air temperature, degrees C. 
% rh: relative humidity, %.
% hu: height of wind measurement, m.
% ht: height of air temperature measurement, m.
% hq: heigh of humidity measurement, m.
% alt: altitude of lake, m.

% Zeng, X., Zhao, M., Dickson, R.E. 1998. Intercomparison of Bulk
% aerodynamic algorithms for the computation of sea surface fluxes using
% TOGA COARE and TAO data. Journal of Climate 11: 2628-2644.

% pre-define arrays
    tstar = nan(length(Uz),1);
    qstar = nan(length(Uz),1);
    u10 = nan(length(Uz),1);
    t10 = nan(length(Uz),1);
    q10 = nan(length(Uz),1);
    wc = nan(length(Uz),1);
    init_wnd = Uz; % initial wind speed
    
% define constants
    const_vonKarman = 0.41; % von Karman constant
    const_gas = 287.1; % gas constant for dry air J kg-1 K-1
    const_SpecificHeatAir = 1005; % Specific heat capacity of air, J kg-1 K-1
    const_Charnock = 0.013; % charnock constant
    const_Gravity = 9.81; % gravitational acceleration, m/s2

% calculate air pressure from altitude
    press = 101325.*(1 - 2.25577e-5.*alt).^5.25588; % Pa
    press = press./100; % mb

% calculate humidity values
    e_s = 6.11.*exp(17.27.*ta./(237.3 + ta)); % saturated vapour pressure at ta, mb
    e_a = rh.*e_s./100; % vapour pressure, mb
    q_z = 0.622.*e_a./press; % specific humidity, kg kg-1 
    e_sat = 6.11.*exp(17.27.*ts./(237.3 + ts)); % saturated vapour pressure at ts, mb
    q_s = 0.622.*e_sat./press; % humidity at saturation, kg kg-1
    
% calculate gas constant for moist air, J kg-1 K-1
    R_a = 287.*(1 + 0.608.*q_z);
    
% calculate latent heat of vaporization, J kg-1 
    xlv = 2.501e6-2370.*ts;    

% calculate air density, kg/m3
    rho_a = 100*press./(R_a.*(ta + 273.16));
    
% kinematic viscosity, m2 s-1
    KinV = (1./rho_a).*(4.94e-8.*ta + 1.7184e-5);    
    
% calculate virtual air temperature, K
    t_virt = (ta + 273.16).*(1 + 0.61.*q_z);     
    
% estimate initial values of u* and zo
    ustar = Uz.*sqrt(0.00104+0.0015./(1+exp((-Uz+12.5)./1.56)));
    zo = (const_Charnock.*ustar.^2./const_Gravity) + (0.11.*KinV./ustar);
    zo_prev = zo.*1.1;
    for i = 1:length(Uz)
        while abs((zo(i) - zo_prev(i)))/abs(zo_prev(i)) > 0.00001
            ustar(i) = const_vonKarman.*Uz(i)/(log(hu./zo(i)));
            dummy = zo(i);
            zo(i)=(const_Charnock.*ustar(i).^2./const_Gravity) + (0.11*KinV(i)./ustar(i));
            zo_prev(i) = dummy;
        end
    end
    
% calculate neutral transfer coefficients
    C_DN = (ustar.^2)./(Uz.^2);
    re = ustar.*zo./KinV;
    zot = zo.*exp(-2.67.*(re).^(0.25) + 2.57);
    zot = real(zot);
    zoq = zot;
    C_HN = const_vonKarman*sqrt(C_DN)./(log(hu./zot)); 
    C_EN = C_HN;
    
% calculate neutral transfer coefficients at 10 m
    C_D10N = (const_vonKarman./log(10./zo)).*(const_vonKarman./log(10./zo)); 
    C_E10N = (const_vonKarman.*const_vonKarman)./...
        (log(10./zo).*log(10./zoq));
    C_H10N = C_E10N;    
    
% calculate neutral latent and sensible heat fluxes, W/m2
    alhN = rho_a.*xlv.*C_EN.*Uz.*(q_s-q_z);
    ashN = rho_a.*const_SpecificHeatAir.*C_HN.*Uz.*(ts-ta); 
    
% calculate initial monin obukhov length scale, m
    obu = (-rho_a.*t_virt.*(ustar.*ustar.*ustar))./...
        (const_vonKarman.*const_Gravity.*(ashN./const_SpecificHeatAir + ...
        0.61.*(ta + 273.16).*alhN./xlv)); 
    
% iteration to compute corrections for atmospheric stability    
    zeta_thres = 15; 
    
    for i = 1:20;          
        
        % calulate roughness lengths
        zo = (0.013.*((ustar.^2)./const_Gravity)) + (0.11.*(KinV./ustar));                        
        re = (ustar.*zo)./KinV;     
        xq = 2.67.*re.^0.25 - 2.57;
        xq(xq < 0) = 0;
        zoq = zo./exp(xq);        
        zot = zoq;    
        
        % define zeta
        zetam = -1.574;
        zetat = -0.465;
        
        % calculate ustar
        zeta = hu./obu;
        zeta(zeta < -zeta_thres) = -zeta_thres;zeta(zeta > zeta_thres) = zeta_thres;
        
        idx = zeta < zetam; % very unstable conditions      
        ustar(idx) = (Uz(idx).*const_vonKarman)./...
            ((log((zetam.*obu(idx))./zo(idx)) - psi_zeng(1,zetam)) + ...
            1.14.*(((-zeta(idx)).^0.333) - ((-zetam).^0.333))); 
        idx = zeta < 0 & zeta >= zetam; % unstable conditions
        ustar(idx) = (Uz(idx).*const_vonKarman)./...
            (log(hu./zo(idx)) - psi_zeng(1,zeta(idx)));
        idx = zeta > 0 & zeta <= 1; % stable conditions
        ustar(idx) = (Uz(idx).*const_vonKarman)./...
            (log(hu./zo(idx)) + 5.*zeta(idx));
        idx = zeta > 1; % very stable conditions
        ustar(idx) = (Uz(idx).*const_vonKarman)./...
            ((log(obu(idx)./zo(idx)) + 5) + (5.*log(zeta(idx)) + zeta(idx) - 1));
        
        % calculate tstar
        zeta = ht./obu;
        zeta(zeta < -zeta_thres) = -zeta_thres;zeta(zeta > zeta_thres) = zeta_thres;
        
        idx = zeta < zetat;% very unstable conditions
        tstar(idx) = (const_vonKarman.*(ta(idx) - ts(idx)))./...
            ((log((zetat.*obu(idx))./zot(idx)) - psi_zeng(2,zetat)) + ...
            0.8.*((-zetat).^-0.333 - ((-zeta(idx))).^-0.333));
        idx = zeta >= zetat & zeta < 0;% unstable conditions
        tstar(idx) = (const_vonKarman.*(ta(idx) - ts(idx)))./...
            (log(ht./zot(idx)) - psi_zeng(2,zeta(idx)));
        idx = zeta > 0 & zeta <= 1; % stable conditions
        tstar(idx) = (const_vonKarman.*(ta(idx) - ts(idx)))./...
            (log(ht./zot(idx)) + 5.*zeta(idx));
        idx = zeta > 1; % very stable conditions
        tstar(idx) = (const_vonKarman.*(ta(idx) - ts(idx)))./...
            ((log(obu(idx)./zot(idx)) + 5) + (5.*log(zeta(idx)) + zeta(idx) - 1));
        
        % calculate qstar
        zeta = hq./obu;
        zeta(zeta < -zeta_thres) = -zeta_thres;zeta(zeta > zeta_thres) = zeta_thres;
        
        idx = zeta < zetat; % very unstable conditions
        qstar(idx) = (const_vonKarman.*(q_z(idx) - q_s(idx)))./...
            ((log((zetat.*obu(idx))./zoq(idx)) - psi_zeng(2,zetat)) + ...
            0.8.*((-zetat).^-0.333 - ((-zeta(idx))).^-0.333));
        idx = zeta >= zetat & zeta < 0; % unstable conditions
        qstar(idx) = (const_vonKarman.*(q_z(idx) - q_s(idx)))./...
            (log(hq./zoq(idx)) - psi_zeng(2,zeta(idx)));
        idx = zeta > 0 & zeta <= 1; % stable conditions
        qstar(idx) = (const_vonKarman.*(q_z(idx) - q_s(idx)))./...
            (log(hq./zoq(idx)) + 5.*zeta(idx));
        idx = zeta > 1; % very stable conditions
        qstar(idx) = (const_vonKarman.*(q_z(idx) - q_s(idx)))./...
            ((log(obu(idx)./zoq(idx)) + 5) + (5.*log(zeta(idx)) + zeta(idx) - 1));
         
        % calculate zeta at 10 m
        zeta = 10./obu;
        zeta(zeta < -zeta_thres) = -zeta_thres;zeta(zeta > zeta_thres) = zeta_thres;
        
        % calculate wind speed at 10 m   
        idx = zeta < zetam; % very unstable conditions 
        u10(idx) = (ustar(idx)./const_vonKarman).*...
            ((log((zetam.*obu(idx))./zo(idx)) - psi_zeng(1,zetam)) + ...
            1.14.*(((-zeta(idx)).^0.333) - ((-zetam).^0.333)));
        idx = zeta < 0 & zeta > zetam;% unstable conditions
        u10(idx) = (ustar(idx)./const_vonKarman).*...
            (log(10./zo(idx)) - psi_zeng(1,zeta(idx)));
        idx = zeta >= 0 & zeta <= 1;% stable conditions    
        u10(idx) = (ustar(idx)./const_vonKarman).*...
            (log(10./zo(idx)) + 5.*zeta(idx));
        idx = zeta > 1; % very stable conditions
        u10(idx) = (ustar(idx)./const_vonKarman).*...
            ((log(obu(idx)./zo(idx)) + 5) + (5.*log(zeta(idx)) + zeta(idx) - 1));
        
        % calcuate air temperature at 10 m
        idx = zeta < zetat;% very unstable conditions 
        t10(idx) = ((tstar(idx)./const_vonKarman).*...
            ((log((zetat.*obu(idx))./zot(idx)) - psi_zeng(2,zetat)) + ...
            0.8.*((-zetat).^-0.333 - ((-zeta(idx))).^-0.333))) + ts(idx);  
        idx = zeta > zetat & zeta < 0;% unstable conditions
        t10(idx) = ((tstar(idx)./const_vonKarman).*...
            log(10./zot(idx)) - psi_zeng(2,zeta(idx))) + ts(idx); 
        idx = zeta >= 0 & zeta <= 1;% stable conditions
        t10(idx) = ((tstar(idx)./const_vonKarman).*...
            (log(10./zot(idx)) + 5.*zeta(idx))) + ts(idx); 
        idx = zeta > 1; % very stable conditions
        t10(idx) = ((tstar(idx)./const_vonKarman).*...
            ((log(obu(idx)./zot(idx)) + 5) + (5.*log(zeta(idx)) + zeta(idx) - 1))) + ts(idx);
        
        % calcuate specific humidity at 10 m
        idx = zeta < zetat; % very unstable conditions
        q10(idx) = ((qstar(idx)./const_vonKarman).*...
            ((log((zetat.*obu(idx))./zoq(idx)) - psi_zeng(2,zetat)) + ...
            0.8.*((-zetat).^-0.333 - ((-zeta(idx))).^-0.333))) + q_s(idx);            
        idx = zeta > zetat & zeta < 0;% unstable conditions
        q10(idx) = ((qstar(idx)./const_vonKarman).*...
            (log(10./zoq(idx)) - psi_zeng(2,zeta(idx)))) + q_s(idx);
        idx = zeta >= 0 & zeta <= 1;% stable conditions
        q10(idx) = ((qstar(idx)./const_vonKarman).*...
            (log(10./zoq(idx)) + 5.*zeta(idx))) + q_s(idx);
        idx = zeta > 1;% very stable conditions
        q10(idx) = ((qstar(idx)./const_vonKarman).*...
            ((log(obu(idx)./zoq(idx)) + 5) + (5.*log(zeta(idx)) + zeta(idx) - 1))) + q_s(idx);

        % calculate transfer coefficients corrected for atmospheric stability
        C_H = (-rho_a.*const_SpecificHeatAir.*ustar.*tstar)./...
            (rho_a.*const_SpecificHeatAir.*Uz.*(ts - ta));
        C_E = C_H;
        C_D = (ustar.*ustar)./(Uz.*Uz);
        
        % calculate tau and sensible and latent heat fluxes
        tau = C_D.*rho_a.*ustar.*ustar;
        ash = rho_a.*const_SpecificHeatAir.*C_H.*Uz.*(ts - ta);           
        alh = rho_a.*xlv.*C_E.*Uz.*(q_s - q_z);
                
        % calculate transfer coefficients at 10m
        C_H10 = ash./(rho_a.*const_SpecificHeatAir.*u10.*(ts - t10));
        C_E10 = alh./(rho_a.*xlv.*u10.*(q_s - q10));
        C_D10 = (ustar.*ustar)./(u10.*u10);
                        
        % calculate new monin obukhov length
        obu = (-rho_a.*t_virt.*(ustar.*ustar.*ustar))./...
            (const_Gravity.*const_vonKarman.*((ash./const_SpecificHeatAir) + ...
            (0.61.*(ta + 273.16).*alh./xlv)));
        
        % alter zeta in stable cases
        zeta = hu./obu;
        idx = zeta >= 1;
        Uz(idx) = max(Uz(idx),0.1);
                        
        % avoid singularity at um = 0 for unstable conditions     
        idx = zeta < 0;
        th = (ta + 273.16).*(1000./press).^(const_gas./1004.67); % potential temperature  
        thvstar = tstar.*(1 + 0.61.*q_z./1000) + 0.61.*th.*qstar; % temperature scaling parameter
        thv = th.*(1 + 0.61.*q_z./1000); % virtual potential temperature    
        wc(idx) = 1.*(-const_Gravity.*ustar(idx).*thvstar(idx)./thv(idx)).^0.333;
        Uz(idx) = sqrt(Uz(idx).*Uz(idx) + wc(idx).*wc(idx));
                
    end
   
% if the measurements are taken at a height of 10 m, keep original input data
    if hu == 10; u10 = init_wnd;end;
    if ht == 10; t10 = ta;end;
       
% convert specific humidity to relative humidity
    if hq == 10; 
        rh10 = rh;
    else
        es = 6.1121.*exp(17.502.*t10./(t10 + 240.97)).*(1.0007+3.46e-6.*press);
        em = q10.*press./(0.378.*q10 + 0.622);
        rh10 = 100.*em./es;
        rh10(rh10 > 100) = 100;
        rh10(rh10 < 0) = 0;
    end
   
% calculate evaporation [mm/day]
    rho_w = 1000*(1-1.9549*0.00001*abs(ts-3.84).^1.68);
    Evap = 86400.*1000.*alh(:)./(rho_w.*xlv);
        
% define output matrix
    mm = [tau(:) alh(:) ash(:) ustar(:) tstar(:) qstar(:) u10(:) t10(:) q10(:),...
        rh10(:) zo(:) zot(:) zoq(:) C_D(:) C_E(:) C_H(:) C_D10(:) C_E10(:) C_H10(:),...
        C_D10N(:) C_E10N(:) C_H10N(:) C_DN(:) C_EN(:) C_HN(:) zeta(:) Evap(:) rho_a(:) q_s(:),...
        ts(:) ta(:) q_z(:) rho_w(:) xlv(:) obu(:)]; 
    mm = real(mm);
end

function psi_zeng = psi_zeng(k,zeta)
    chik = (1 - 16.*zeta).^0.25;
    if k == 1
        psi_zeng = 2.*log((1 + chik).*0.5) + log((1+chik.*chik).*0.5) - ...
            2.*atan(chik) + (pi./2);
    else
        psi_zeng = 2.*log((1 + chik.*chik).*0.5);
    end
end