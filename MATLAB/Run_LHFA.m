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
    
    % remove nand
    idx = isnan(wtr);
    wtr(idx) = [];
    wtrD(idx) = []; 
    fprintf('...completed\n\n') ;
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
    
    fprintf('...completed\n\n') ;
end   

if TT.openSW
    fprintf(['Reading ' LakeName '.sw file'])
    swFileName  = [Folder '/' LakeName '.sw'];
    oper = fopen(swFileName);
    if eq(oper,-1)
        error([LakeName '.sw file not found']);
    end
    fclose all;
    [swD,sw] = gFileOpen(swFileName);
    sw(sw < 0) = 0;
    
    % remove nans
    idx = isnan(sw);
    sw(idx) = [];
    swD(idx) = []; 
    fprintf('...completed\n\n') ;
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
    fprintf('...completed\n\n') ;
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
end

% make sure dates match
if ~TT.openSW
    if TT.openWtr && TT.openWnd && TT.openRH && TT.openAirT 
        idx = intersect(intersect(intersect(wtrD,wndD),rhD),airTD); 
        dates = idx;
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
    end
end

%*** find samplerate of raw data 
if length(dates) < smplTs
    tLen = length(dates);
else
    tLen = smplTs;
end
steps = NaN(1,tLen-1);
for i = 1:tLen-1
    steps(i) = dates(i+1)-dates(i);
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
    fprintf('Down sampling data');
    if TT.openWtr
        [DS_wtrD,DS_wtr] = DownSample_TS(wtrD,outRs,wtr);
        wtrD = DS_wtrD;
        wtr = DS_wtr;
        varL = length(wtrD);
        dates = wtrD;
        clear DS_wtrD DS_wtr
    end
    if TT.openAirT
        [DS_airTD,DS_airT] = DownSample_TS(airTD,outRs,airT);
        airTD = DS_airTD;
        airT = DS_airT;
        varL = length(airTD);
        dates = airTD;
        clear DS_airTD DS_airT
    end
    if TT.openRH
        [DS_rhD,DS_rh] = DownSample_TS(rhD,outRs,rh);
        rhD = DS_rhD;
        rh = DS_rh;
        varL = length(rhD);
        dates = rhD;
        clear DS_rhD DS_rh
    end
    if TT.openSW
        [DS_swD,DS_sw] = DownSample_TS(swD,outRs,sw);
        swD = DS_swD;
        sw = DS_sw;
        varL = length(swD);
        dates = swD;
        clear DS_swD DS_sw
    end
    if TT.openWnd
        [DS_wndD,DS_wnd] = DownSample_TS(wndD,outRs,wnd);
        wndD = DS_wndD;
        wnd = DS_wnd;
        varL = length(wndD);
        dates = wndD;
        clear DS_swD DS_sw
    end
    fprintf('...completed\n\n');
end

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
end
    
%****-----varL is the length of output files as of here-------*****
% water temperature
if TT.wrt_wTemp
    writeTable.wTemp = wtr;
end

% calculate surface fluxes
if TT.senslatYes
    mm = sens_latent(wtr,wnd,airT,rh,wndH,htH,hqH,alt);
end

% monin-obukhow length scale
if TT.wrt_obu
    writeTable.obu = mm(:,35);
end

% momentum flux
if TT.wrt_tau
    writeTable.tau = mm(:,1);
end

% sensible heat flux
if TT.wrt_Qh
    writeTable.Qh = mm(:,3);
end

% latent heat flux
if TT.wrt_Qe
    writeTable.Qe = mm(:,2);
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
    press = 101325.*(1 - 2.25577e-5.*alt).^5.25588; % Pa
    press = press./100; % mb
    [~,~,Qlnet] = calc_lwnet(wtrD,lat,press,airT,rh,sw,wtr);
    
    % ensure night time long wave fluxes are equal to the average day time fluxes
    % cannot detect clooud cover during night
    dateV = datevec(dates);
    [~,~,b] = unique(dateV(:,1:3),'rows');
    daily_av = accumarray(b,Qlnet,[],@nanmean);
    Q2 = Qlnet;
    for ii = 1:length(unique(b));
        Q2(b == ii) = daily_av(ii);
    end
    % then replace nan with dat2
    Qlnet(isnan(Qlnet)) = Q2(isnan(Qlnet));
    
    writeTable.Qlnet = -Qlnet;
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
    
    writeTable.Qlin = -lw;
end

% outgoing long wave heat flux
if TT.wrt_Qlout
    Tk = wtr + 273.13;
    emiss = 0.972;
    S_B = 5.67E-8;
    LWo = S_B*emiss*Tk.^4;
    writeTable.Qlout = -LWo;
end

% reflected short wave radiaiton
if TT.wrt_Qsr
    sw_alb = sw_albedo(dates,lat);
    Qsr = sw.*sw_alb; % reflected short wave radiation
    writeTable.Qsr = Qsr;
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
        wrt([delimO writeNames{i}]);
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
