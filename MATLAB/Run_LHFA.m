function Run_LHFA(LakeName,Folder,skipLoad)
%----Author: Jordan S Read 2009 ----
%----version 3.3.1 2013-05-10 ------

if nargin < 3
    skipLoad = false;
end
clc; close all
done = false;
if ~skipLoad
    build_config(LakeName,Folder);
    while ~done
        pause(.4)
    end
    pause(0.1);
end

% -- variables --
matSec = 86400;         % number of seconds in a day
smplTs = 50;            % samplerate test length
dateTl = 0.00001;       % tolerance for equal day (day fraction)
fprintf(['Reading ' LakeName '.hfx file']);
[outPuts,outRs,wndH,hqH,htH,lat,alt,...
            wndMx,wndMn,plotYes,writeYes] = ...
    OpenCfg( LakeName,Folder );
% -- end variables --

if ~any([plotYes writeYes])    
    error(['User must specify either to write results,'...
        ' plot results, or both ....... (see ' LakeName '.hfx)'])
end

fprintf('...completed\n\n') ;

fprintf('****Building program structure****\n');
pltMods = [];
if plotYes
    fprintf(['Checking for ' LakeName '.plt file*'])
    pltFileName  = [Folder '/' LakeName '.plt'];
    oper = fopen(pltFileName);
    if eq(oper,-1)
        fprintf('...not found [*is optional]\n\n')
    else
        
        pltMods = pltFileOpen(pltFileName);
        fprintf('...completed\n');
    end
end
[TT,outputOptions,writeTable,plotTable,~,dateOutput,~,delimO] = ...
    OutputConstructor(outPuts,pltMods);
fNms = fieldnames(TT);
for j = 1:length(fNms)
    if TT.(char(fNms{j}))
       fprintf([char(fNms{j}) '\n'])
    end
end
fprintf('****completed****\n\n') ;

if TT.openWtr
    fprintf(['Reading ' LakeName '.wtr file'])
    wtrFileName  = [Folder '/' LakeName '.wtr'];
    oper = fopen(wtrFileName);
    if eq(oper,-1)
        error([LakeName '.wtr file not found']);
    end
    fclose all;
    [wtrD,wtr,heads] = gFileOpen(wtrFileName);
    headers = textscan(heads,'%s','Delimiter','\t');
    headers = headers{1}(2:end); % get rid of dateTime
    depths = NaN(1,length(headers));
    for d =1:length(depths)
        txt = headers{d};
        depths(d) = str2double(txt(5:end));
    end
    [mnDep,mnIdx] = min(depths);
    if length(depths) > 1;
        wtr = wtr(:,mnIdx);
    end  
    if mnDep>0
        disp([' ...' num2str(mnDep) ...
            'm is the shallowest depth in the .wtr file'...
            ' which will be used to represent surface water temperatures'])
        pause(1.5);
    end
    
    % remove nans
    idx = isnan(wtr);
    wtr(idx) = [];
    wtrD(idx) = []; 
    fprintf('...completed\n\n') ;
    
    %*** find samplerate of raw data 
    if length(wtrD) < smplTs
        tLen = length(wtrD);
    else
        tLen = smplTs;
    end
    steps = NaN(1,tLen-1);
    for i = 1:tLen-1
        steps(i) = wtrD(i+1)-wtrD(i);
    end    
    if eq(min(steps),0)
        matRs = mean(steps)*matSec;
    else
        matRs = min(steps)*matSec; %current sample rate of raw data in seconds
    end
    clear vals ind numMx numCont steps tLen
    if (outRs - matRs)>dateTl
        TT.dwnSmple = true;   %down sample if necessary
    end

    % *** down sampling ***
    if TT.dwnSmple
        fprintf(['Down sampling ' LakeName '.wtr data']);
        [DS_wtrD,DS_wtr] = DownSample_TS(wtrD,outRs,wtr);
        wtrD = DS_wtrD;
        wtr = DS_wtr;
        varL = length(wtrD);
        clear DS_wtrD DS_wtr
        fprintf('...completed\n\n');
    end
end
           
if TT.openWnd
    fprintf(['Reading ' LakeName '.wnd file'])
    wndFileName  = [Folder '/' LakeName '.wnd'];
    oper = fopen(wndFileName);
    if eq(oper,-1)
        error([LakeName '.wnd file not found']);
    end
    fclose all;
    [wndD,wnd] = gFileOpen(wndFileName);
    wnd(wnd < wndMn) = wndMn;
    wnd(wnd > wndMx) = wndMx;    
    
    % remove nans
    idx = isnan(wnd);
    wnd(idx) = [];
    wndD(idx) = []; 
    fprintf('...completed\n\n');
    
    %*** find samplerate of raw data 
    if length(wndD) < smplTs
        tLen = length(wndD);
    else
        tLen = smplTs;
    end
    steps = NaN(1,tLen-1);
    for i = 1:tLen-1
        steps(i) = wndD(i+1)-wndD(i);
    end    
    if eq(min(steps),0)
        matRs = mean(steps)*matSec;
    else
        matRs = min(steps)*matSec; %current sample rate of raw data in seconds
    end
    clear vals ind numMx numCont steps tLen
    if (outRs - matRs)>dateTl
        TT.dwnSmple = true;   %down sample if necessary
    end

    % *** down sampling ***
    if TT.dwnSmple
        fprintf(['Down sampling ' LakeName '.wnd data']);
        [DS_wndD,DS_wnd] = DownSample_TS(wndD,outRs,wnd);
        wndD = DS_wndD;
        wnd = DS_wnd;
        varL = length(wndD);
        clear DS_wndD DS_wnd
        fprintf('...completed\n\n');
    end
    dates = wndD;
end   

if TT.openSW
    if exist([Folder '/' LakeName '.sw']) > 0;
        fprintf(['Reading ' LakeName '.sw file'])
        swFileName  = [Folder '/' LakeName '.sw'];
        fclose all;
        [swD,sw] = gFileOpen(swFileName);
        sw(sw < 0) = 0;

        % remove nans
        idx = isnan(sw);
        sw(idx) = [];
        swD(idx) = []; 
        fprintf('...completed\n\n');
        
        %*** find samplerate of raw data 
        if length(swD) < smplTs
            tLen = length(swD);
        else
            tLen = smplTs;
        end
        steps = NaN(1,tLen-1);
        for i = 1:tLen-1
            steps(i) = swD(i+1)-swD(i);
        end    
        if eq(min(steps),0)
            matRs = mean(steps)*matSec;
        else
            matRs = min(steps)*matSec; %current sample rate of raw data in seconds
        end
        clear vals ind numMx numCont steps tLen
        if (outRs - matRs)>dateTl
            TT.dwnSmple = true;   %down sample if necessary
        end

        % *** down sampling ***
        if TT.dwnSmple
            fprintf(['Down sampling ' LakeName '.sw data']);
            [DS_swD,DS_sw] = DownSample_TS(swD,outRs,sw);
            swD = DS_swD;
            sw = DS_sw;
            varL = length(swD);
            clear DS_swD DS_sw
            fprintf('...completed\n\n');
        end
        
    else if exist([Folder '/' LakeName '.par']) > 0;
        fprintf(['Reading ' LakeName '.par file'])
        swFileName  = [Folder '/' LakeName '.par'];
        fclose all;
        [swD,sw] = gFileOpen(swFileName);
        sw(sw < 0) = 0;
        % remove nans
        idx = isnan(sw);
        sw(idx) = [];
        swD(idx) = []; 
        parMult = 0.4957;
        sw = sw.*parMult;
        fprintf('...completed\n\n');
                
        %*** find samplerate of raw data 
        if length(swD) < smplTs
            tLen = length(swD);
        else
            tLen = smplTs;
        end
        steps = NaN(1,tLen-1);
        for i = 1:tLen-1
            steps(i) = swD(i+1)-swD(i);
        end    
        if eq(min(steps),0)
            matRs = mean(steps)*matSec;
        else
            matRs = min(steps)*matSec; %current sample rate of raw data in seconds
        end
        clear vals ind numMx numCont steps tLen
        if (outRs - matRs)>dateTl
            TT.dwnSmple = true;   %down sample if necessary
        end

        % *** down sampling ***
        if TT.dwnSmple
            fprintf(['Down sampling ' LakeName '.sw data']);
            [DS_swD,DS_sw] = DownSample_TS(swD,outRs,sw);
            swD = DS_swD;
            sw = DS_sw;
            varL = length(swD);
            clear DS_swD DS_sw
            fprintf('...completed\n\n');
        end
        
        else
            error([LakeName '.sw nor .par file not found']);
        end
    end
    dates = swD;
end  

if TT.openAirT
    fprintf(['Reading ' LakeName '.airT file'])
    airTFileName  = [Folder '/' LakeName '.airT'];
    oper = fopen(airTFileName);
    if eq(oper,-1)
        error([LakeName '.airT file not found']);
    end
    fclose all;
    [airTD,airT] = gFileOpen(airTFileName);
    
    % remove nans
    idx = isnan(airT);
    airT(idx) = [];
    airTD(idx) = []; 
    fprintf('...completed\n\n');
    
    %*** find samplerate of raw data 
    if length(airTD) < smplTs
        tLen = length(airTD);
    else
        tLen = smplTs;
    end
    steps = NaN(1,tLen-1);
    for i = 1:tLen-1
        steps(i) = airTD(i+1)-airTD(i);
    end
    if eq(min(steps),0)
        matRs = mean(steps)*matSec;
    else
        matRs = min(steps)*matSec; %current sample rate of raw data in seconds
    end
    clear vals ind numMx numCont steps tLen
    if (outRs - matRs)>dateTl
        TT.dwnSmple = true;   %down sample if necessary
    end
    
    % *** down sampling ***
    if TT.dwnSmple
        fprintf(['Down sampling ' LakeName '.airT data']);
        [DS_airTD,DS_airT] = DownSample_TS(airTD,outRs,airT);
        airTD = DS_airTD;
        airT = DS_airT;
        varL = length(airTD);
        clear DS_airTD DS_airT
        fprintf('...completed\n\n');
    end
end   

if TT.openRH
    fprintf(['Reading ' LakeName '.rh file'])
    rhFileName  = [Folder '/' LakeName '.rh'];
    oper = fopen(rhFileName);
    if eq(oper,-1)
        error([LakeName '.rh file not found']);
    end
    fclose all;
    [rhD,rh] = gFileOpen(rhFileName);
    rh(rh < 0) = 0;
    rh(rh > 100) = 100;  
    
    % remove nans
    idx = isnan(rh);
    rh(idx) = [];
    rhD(idx) = []; 
    fprintf('...completed\n\n') ;
    
    %*** find samplerate of raw data 
    if length(rhD) < smplTs
        tLen = length(rhD);
    else
        tLen = smplTs;
    end
    steps = NaN(1,tLen-1);
    for i = 1:tLen-1
        steps(i) = rhD(i+1)-rhD(i);
    end
    if eq(min(steps),0)
        matRs = mean(steps)*matSec;
    else
        matRs = min(steps)*matSec; %current sample rate of raw data in seconds
    end
    clear vals ind numMx numCont steps tLen
    if (outRs - matRs)>dateTl
        TT.dwnSmple = true;   %down sample if necessary
    end
    
    % *** down sampling ***
    if TT.dwnSmple
        fprintf(['Down sampling ' LakeName '.rh data']);
        [DS_rhD,DS_rh] = DownSample_TS(rhD,outRs,rh);
        rhD = DS_rhD;
        rh = DS_rh;
        varL = length(rhD);
        clear DS_rhD DS_rh
        fprintf('...completed\n\n');
    end
end

% make sure dates match
if ~TT.openSW
    if TT.openWtr && TT.openWnd && TT.openRH && TT.openAirT 
        idx = intersect(intersect(intersect(wtrD,wndD),rhD),airTD); 
        dates = idx;
        varL = length(dates);
        wtrD = idx;
        wndD = idx;
        rhD = idx;
        airTD = idx;
        wtr = wtr(ismember(wtrD,idx));
        wnd = wnd(ismember(wndD,idx));
        rh = rh(ismember(rhD,idx));
        airT = airT(ismember(airTD,idx)); 
    end
end

if ~TT.openWnd
    if TT.openRH && TT.openAirT && TT.openWtr && TT.openSW 
        idx = intersect(intersect(intersect(rhD,airTD),wtrD),swD);    
        dates = idx;
        varL = length(dates);
        wtrD = idx;
        swD = idx;
        rhD = idx;
        airTD = idx;
        rh = rh(ismember(rhD,idx));
        airT = airT(ismember(airTD,idx));
        wtr = wtr(ismember(wtrD,idx));
        sw = sw(ismember(swD,idx)); 
    end
end

if TT.openRH && TT.openAirT && TT.openWtr && TT.openSW && TT.openWnd
    idx = intersect(intersect(intersect(intersect(rhD,airTD),wtrD),swD),wndD);    
    dates = idx;
    varL = length(dates);
    wtrD = idx;
    swD = idx;
    rhD = idx;
    airTD = idx;
    wndD = idx;
    rh = rh(ismember(rhD,idx));
    airT = airT(ismember(airTD,idx));
    wtr = wtr(ismember(wtrD,idx));
    sw = sw(ismember(swD,idx)); 
    wnd = wnd(ismember(wndD,idx)); 
end

if ~TT.openRH && ~TT.openAirT && ~TT.openSW && ~TT.openWnd
    if TT.openWtr   
        dates = wtrD;
        varL = length(dates);
    end
end

% look for long-wave radiation data
if TT.openLWnet;
    if exist([Folder '/' LakeName '.lwnet']) > 0;
        fprintf(['Reading ' LakeName '.lwnet file'])
        lwnetFileName = [Folder '/' LakeName '.lwnet'];
        fclose all;
        [lwnetD,lwnet] = gFileOpen(lwnetFileName);
        idx = isnan(lwnet);
        lwnet(idx) = [];
        lwnetD(idx) = [];      
        fprintf('...completed\n\n');
        
        %*** find samplerate of raw data 
        if length(lwnetD) < smplTs
            tLen = length(lwnetD);
        else
            tLen = smplTs;
        end
        steps = NaN(1,tLen-1);
        for i = 1:tLen-1
            steps(i) = lwnetD(i+1)-lwnetD(i);
        end
        if eq(min(steps),0)
            matRs = mean(steps)*matSec;
        else
            matRs = min(steps)*matSec; %current sample rate of raw data in seconds
        end
        clear vals ind numMx numCont steps tLen
        if (outRs - matRs)>dateTl
            TT.dwnSmple = true;   %down sample if necessary
        end

        % *** down sampling ***
        if TT.dwnSmple
            fprintf(['Down sampling ' LakeName '.lwnet data']);
            [DS_lwnetD,DS_lwnet] = DownSample_TS(lwnetD,outRs,lwnet);
            lwnetD = DS_lwnetD;
            lwnet = DS_lwnet;
            varL = length(lwnetD);
            clear DS_lwnetD DS_lwnet
            fprintf('...completed\n\n');
        end
        
    else if exist([Folder '/' LakeName '.lw']) > 0;
            fprintf([LakeName '.lwnet files not found, looking for .lw file...\n\n'])
            fprintf(['...Reading ' LakeName '.lw file']);
            lwFileName  = [Folder '/' LakeName '.lw'];
            fclose all;
            [lwD,lw] = gFileOpen(lwFileName);
            idx = isnan(lw);
            lw(idx) = [];
            lwD(idx) = [];      
            fprintf('...completed\n\n');

            %*** find samplerate of raw data 
            if length(lwD) < smplTs
                tLen = length(lwD);
            else
                tLen = smplTs;
            end
            steps = NaN(1,tLen-1);
            for i = 1:tLen-1
                steps(i) = lwD(i+1)-lwD(i);
            end
            if eq(min(steps),0)
                matRs = mean(steps)*matSec;
            else
                matRs = min(steps)*matSec; %current sample rate of raw data in seconds
            end
            clear vals ind numMx numCont steps tLen
            if (outRs - matRs)>dateTl
                TT.dwnSmple = true;   %down sample if necessary
            end

            % *** down sampling ***
            if TT.dwnSmple
                fprintf(['Down sampling ' LakeName '.lw data']);
                [DS_lwD,DS_lw] = DownSample_TS(lwD,outRs,lw);
                lwD = DS_lwD;
                lw = DS_lw;
                varL = length(lwD);
                clear DS_lwD DS_lw
                fprintf('...completed\n\n');
            end
                    
            % find when wtr and lw dates intersect
            idx = intersect(wtrD,lwD);
            lw = lw(ismember(lwD,idx));
            wtr = wtr(ismember(wtrD,idx));
            wtrD = idx;
            
            Tk = wtr + 273.13; % .wtr already called at this point
            emiss = 0.972;
            S_B = 5.67E-8;
            LWo = S_B*emiss*Tk.^4;
                        
            % define lwnet
            lwnet = lw - LWo;
            lwnet = lwnet;
            lwnetD = idx;
        else
            fprintf(['...' LakeName '.lwnet and .lw file missing, using .airT, .rh, .sw, .wtr instead...\n\n'])
            
            press = 101325.*(1 - 2.25577e-5.*alt).^5.25588; % Pa
            press = press./100; % mb
            [~,~,Qlnet] = calc_lwnet(dates,lat,press,airT,rh,sw,wtr);
            lwnet = Qlnet;
            lwnetD = dates;
        end
    end     
end

% re-adjust dates depending on lw and lwnet files
if TT.openRH && TT.openAirT && TT.openWtr && TT.openSW && TT.openLWnet
    idx = intersect(intersect(intersect(intersect(rhD,airTD),wtrD),swD),lwnetD);    
    dates = idx;
    varL = length(dates);
    wtrD = idx;
    swD = idx;
    rhD = idx;
    airTD = idx;
    lwnetD = idx;
    rh = rh(ismember(rhD,idx));
    airT = airT(ismember(airTD,idx));
    wtr = wtr(ismember(wtrD,idx));
    sw = sw(ismember(swD,idx)); 
    lwnet = lwnet(ismember(lwnetD,idx)); 
end

if TT.openRH && TT.openAirT && TT.openWtr && TT.openSW && TT.openLWnet && TT.openWnd
    idx = intersect(intersect(intersect(intersect(intersect(rhD,airTD),wtrD),swD),lwnetD),wndD);    
    dates = idx;
    varL = length(dates);
    wtrD = idx;
    swD = idx;
    rhD = idx;
    airTD = idx;
    lwnetD = idx;
    wndD = idx;
    rh = rh(ismember(rhD,idx));
    airT = airT(ismember(airTD,idx));
    wtr = wtr(ismember(wtrD,idx));
    sw = sw(ismember(swD,idx)); 
    lwnet = lwnet(ismember(lwnetD,idx)); 
    wnd = wnd(ismember(wndD,idx)); 
end

% if data isn't downsampled, define times and variable length
if ~TT.dwnSmple
    if TT.openWtr
        dates = wtrD;
        varL = length(dates);
    end
    if TT.openAirT
        dates = airTD;
        varL = length(dates);
    end
    if TT.openRH
        dates = rhD;
        varL = length(dates);
    end
    if TT.openSW
        dates = swD;
        varL = length(dates);
    end
    if TT.openWnd
        dates = wndD;
        varL = length(dates);
    end
    if TT.openLWnet
        dates = lwnetD;
        varL = length(dates);
    end
end
    
%****-----varL is the length of output files as of here-------*****
% water temperature
if TT.wrt_wTemp
    writeTable.wTemp = wtr;
end

% calculate surface fluxes
if TT.senslatYes || TT.QtotYes
    mm = sens_latent(wtr,wnd,airT,rh,wndH,htH,hqH,alt,lat);
end

% atmospheric stability
if TT.wrt_obu
    zL1 = wndH./mm(:,35);
    zL1(zL1 > 15) = 15;
    zL1(zL1 < -15) = -15;
    writeTable.obu = zL1;
end

% momentum flux
if TT.wrt_tau
    writeTable.tau = mm(:,1);
end

% sensible heat flux
if TT.wrt_Qh || TT.QtotYes
    Qh = mm(:,3);
    if TT.wrt_Qh
        writeTable.Qh = Qh;
    end
end

% latent heat flux
if TT.wrt_Qe || TT.QtotYes
    Qe = mm(:,2);
    if TT.wrt_Qe
        writeTable.Qe = Qe;
    end
end

% air shear velocity 
if TT.wrt_uSt_a
    writeTable.uSt_a = mm(:,4);
end

% air shear velocity (neutral) 
if TT.wrt_uSt_aN
    mm2 = neutral_transfer_coeff(wnd,wndH);
    writeTable.uSt_aN = mm2(:,1);
end

% wind speed at 10 m
if TT.wrt_u10
    writeTable.u10 = mm(:,7);
end

% wind speed at 10 m (neutral)
if TT.wrt_u10N
    mm2 = neutral_transfer_coeff(wnd,wndH);
    writeTable.u10N = mm2(:,2);
end

% air temperature at 10 m
if TT.wrt_t10
    writeTable.t10 = mm(:,8);
end

% relative humidity at 10 m
if TT.wrt_rh10
    writeTable.rh10 = mm(:,10);
end

% transfer coefficient for momentum
if TT.wrt_C_D
    writeTable.C_D = mm(:,14);
end

% transfer coefficient for heat
if TT.wrt_C_E
    writeTable.C_E = mm(:,15);
end

% transfer coefficient for humidity
if TT.wrt_C_H
    writeTable.C_H = mm(:,16);
end

% transfer coefficient for momentum at 10 m
if TT.wrt_C_D10
    writeTable.C_D10 = mm(:,17);
end

% transfer coefficient for heat at 10 m
if TT.wrt_C_E10
    writeTable.C_E10 = mm(:,18);
end

% transfer coefficient for humidity at 10 m
if TT.wrt_C_H10
    writeTable.C_H10 = mm(:,19);
end

% neutral transfer coefficient for momentum
if TT.wrt_C_D10N
    mm2 = neutral_transfer_coeff(wnd,wndH);
    writeTable.C_D10N = mm2(:,3);
end

% neutral drag coefficient for heat
if TT.wrt_C_E10N
    mm2 = neutral_transfer_coeff(wnd,wndH);
    writeTable.C_E10N = mm2(:,4);
end

% neutral drag coefficient for humidity
if TT.wrt_C_H10N
    mm2 = neutral_transfer_coeff(wnd,wndH);
    writeTable.C_H10N = mm2(:,5);
end

% neutral transfer coefficient for momentum
if TT.wrt_C_DN
    mm2 = neutral_transfer_coeff(wnd,wndH);
    writeTable.C_DN = mm2(:,6);
end

% neutral drag coefficient for heat
if TT.wrt_C_EN
    mm2 = neutral_transfer_coeff(wnd,wndH);
    writeTable.C_EN = mm2(:,7);
end

% neutral drag coefficient for humidity
if TT.wrt_C_HN
    mm2 = neutral_transfer_coeff(wnd,wndH);
    writeTable.C_HN = mm2(:,8);
end

% evaporation
if TT.wrt_Evap
    writeTable.Evap = mm(:,27);
end
 
% net long wave heat flux
if TT.wrt_Qlnet
    writeTable.Qlnet = lwnet;
end

% incoming long wave heat flux
if TT.wrt_Qlin
    press = 101325.*(1 - 2.25577e-5.*alt).^5.25588; % Pa
    press = press./100; % mb
    [lw,~,~] = calc_lwnet(wtrD,lat,press,airT,rh,sw,wtr);
    
    % ensure night time long wave fluxes are equal to the average day time fluxes
    dateV = datevec(dates);
    [~,~,b] = unique(dateV(:,1:3),'rows');
    daily_av = accumarray(b,lw,[],@nanmean);
    Q2 = lw;
    for ii = 1:length(unique(b));
        Q2(b == ii) = daily_av(ii);
    end
    % then replace nan with dat2
    lw(isnan(lw)) = Q2(isnan(lw));
    
    writeTable.Qlin = lw;
end

% outgoing long wave heat flux
if TT.wrt_Qlout
    Tk = wtr + 273.13;
    emiss = 0.972;
    S_B = 5.67E-8;
    LWo = S_B*emiss*Tk.^4;
    writeTable.Qlout = LWo;
end

% reflected short wave radiaiton
if TT.wrt_Qsr || TT.QtotYes
    sw_alb = sw_albedo(dates,lat);
    Qsr = sw.*sw_alb; % reflected short wave radiation
    
    if TT.wrt_Qsr
        writeTable.Qsr = Qsr;
    end
end

% total surface heat flux
if TT.wrt_Qtot
    Qtot = sw - Qsr - Qe - Qh + lwnet;
    writeTable.Qtot = Qtot;
end

% build plot array
if plotYes
    fprintf('Plotting results');
    plotLA_results(writeTable,plotTable,dates,LakeName,Folder)
    fprintf('...completed\n\n');
end

% build file array
writeNames = {};
cnt = 1;
for k = 1:length(outputOptions)
    if ~islogical(writeTable.(char(outputOptions{k})))
        writeNames{cnt} = outputOptions{k};
        cnt = cnt+1;
    end
end

% write to file
if writeYes
    fprintf('Writing results to file');
end

if writeYes && ~isempty(writeNames)
    outputFile = [Folder '/' LakeName '_results.txt'];
    outFile = fopen(outputFile,'w');
    if eq(outFile,-1)
        error([Folder '/' LakeName '_results.csv file in use, please close']);
    end
    wrt = @(writer)fprintf(outFile,writer); % build a subfunction that writes 
                                    % the contents of the input "writer" 
                                    % to the file everytime wrt is called
    wrt('DateTime');
    for i = 1:cnt-1
        wrt([delimO writeNames{i} plotTable.(writeNames{i}).YLabel]);
    end
    wrt('\r\n');
    for j = 1:varL
        wrt(datestr(dates(j),dateOutput)); %change 'dateOutput' 
                                    % in the 'OutputConstructor.m' file
        for i = 1:length(writeNames)
            wrt([delimO num2str(writeTable.(char(writeNames{i}))(j))]);
        end
        wrt('\r\n');
    end
    fclose all;
end
if writeYes
    fprintf('...completed\n\n');
end
disp('Lake Heat Flux Analyzer is complete')
%profile off
%profile viewer
end
