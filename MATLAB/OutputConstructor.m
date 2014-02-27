function [truthTable,outputOptions,writeTable,plotTable,...
    dateInput,dateOutput,delimI,delimO] ...
    = OutputConstructor(outputNames,pltMods)

outputOptions = {'tau','Qh','Qe','C_DN','C_EN',...
     'C_HN','C_D10N','C_E10N','C_H10N','C_D','C_E','C_H',...
     'C_D10','C_E10','C_H10','u10','u10N','t10','rh10','Qlnet','Qlin','Qlout',...
     'uSt_a','uSt_aN','Evap','wTemp','Qsr','obu','Qtot'};
  
% writable outputs
programOptions= {'openWtr','openWnd','openWnd','openRH',...
    'openAirT','openSW','openLWnet','dwnSmple','senslatYes','QtotYes'};
% program flow specifiers

% Error Check
if ne(length(outputNames),sum(ismember(outputNames,outputOptions)))
    error(['output "' ...
        outputNames{~ismember(outputNames,outputOptions)} ...
        '" not recognized'])
end
% Other defaults
dateInput  = 'yyyy-mm-dd HH:MM';
delimI = '\t';

dateOutput  = 'yyyy-mm-dd HH:MM';
delimO = '\t';

% Figure defaults
isStringMod = {'figUnits','figType','fontName','figRes'};
figUnits    = 'inches';
figWidth    = 6;   % relative to fig_units
figHeight   = 3;   % relative to fig_units
leftMargin  = .75; % relative to fig_units
rightMargin = .1;  % relative to fig_units
topMargin   = .4;  % relative to fig_units
botMargin   = .4;  % relative to fig_units
figType     = 'png';
figRes      = '150'; % dots per inch (not relative to units?)
fontName    = 'Arial';
fontSize    = 12;
heatMapMin  = 0;
heatMapMax  = 30;


plt = struct('figUnits',figUnits,'figWidth',figWidth,'figHeight',figHeight,...
    'leftMargin',leftMargin,'rightMargin',rightMargin,'topMargin',topMargin,...
    'botMargin',botMargin,'figType',figType,'figRes',figRes,...
    'fontName',fontName,'fontSize',fontSize,'heatMapMin',heatMapMin,...
    'heatMapMax',heatMapMax);

if ~isempty(pltMods)
    % use plot mods to modify plotting defaults
    fN = fieldnames(pltMods);
    for n = 1:length(fN)
        if ~any(strcmp(fN{n},isStringMod))
            if ne(str2double(pltMods.(fN{n})),plt.(fN{n}))
                fprintf(['>>User plot modification; replacing ' fN{n} '=' ...
                    num2str(plt.(fN{n})) ' with ' pltMods.(fN{n})  '<<<\n']);
                plt.(fN{n}) = str2double(pltMods.(fN{n}));
            end
        else
            if ~strcmp(pltMods.(fN{n}),plt.(fN{n}))
                fprintf(['>>User plot modification; replacing ' fN{n} '=' ...
                    plt.(fN{n}) ' with ' pltMods.(fN{n})  '<<<\n']);
                plt.(fN{n}) = pltMods.(fN{n});
            end
        end
    end
end


fig_Defaults = struct('Units',plt.figUnits,'Color','w',...
    'PaperUnits',plt.figUnits,...
    'PaperPosition',[0 0 plt.figWidth plt.figHeight],...
    'Position',[1 1 plt.figWidth plt.figHeight],...
    'PaperPositionMode','manual',...
    'PaperSize',[plt.figWidth plt.figHeight]);

print_Defaults = struct('format',['-d' plt.figType],'res',['-r' plt.figRes],...
    'toClose',true);

position = [plt.leftMargin/plt.figWidth plt.botMargin/plt.figHeight ...
    (plt.figWidth-plt.leftMargin-plt.rightMargin)/plt.figWidth ...
    (plt.figHeight-plt.topMargin-plt.botMargin)/plt.figHeight];

axes_Defaults = struct('FontName',plt.fontName,'FontSize',plt.fontSize,...
    'Layer','top','Position',position,...
    'Box','on',...
    'YLabel','','Title','','YDir','normal','YScale','linear','CLim',[0 1]);
% build Structures
plotTable = struct('FigD',fig_Defaults,'PrintD',print_Defaults);
for j = 1:length(programOptions)
    truthTable.(char(programOptions{j})) = false;
end
for j = 1:length(outputOptions)
    truthTable.(['wrt_' char(outputOptions{j})]) = false;
    writeTable.(char(outputOptions{j})) = false;
    plotTable.(char(outputOptions{j})) = axes_Defaults;
    
end
for j = 1:length(outputOptions)
    outputConstruct.(char(outputOptions{j})) = truthTable;
end

% Water Temperature
name = 'wTemp';
RunNeed.(char(name)) = {'openWtr','wrt_wTemp'};
RunAxes.(char(name)) = struct('YLabel',' (^{o} C)',...
    'Title','Surface water temperature');

% Momentum flux
name = 'tau';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_tau',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel',' (N m^{-2})',...
    'Title','Surface flux of momentum');

% Sensible heat flux
name = 'Qh';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_Qh',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel',' (W m^{-2})',...
    'Title','Sensible heat flux');

% Latent heat flux
name = 'Qe';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_Qe',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel',' (W m^{-2})',...
    'Title','Latent heat flux');

% Transfer coefficient for momentum (neutral)
name = 'C_DN';
RunNeed.(char(name)) = {'openWnd','wrt_C_DN'};
RunAxes.(char(name)) = struct('YLabel','',...
    'Title','Transfer coefficient for momentum (neutral)');

% Transfer coefficient for heat (neutral)
name = 'C_EN';
RunNeed.(char(name)) = {'openWnd','wrt_C_EN'};
RunAxes.(char(name)) = struct('YLabel','',...
    'Title','Transfer coefficient for heat (neutral)');

% Transfer coefficient for humidity (neutral)
name = 'C_HN';
RunNeed.(char(name)) = {'openWnd','wrt_C_HN'};
RunAxes.(char(name)) = struct('YLabel','',...
    'Title','Transfer coefficient for humidity (neutral)');

% Transfer coefficient for momentum at 10 m (neutral)
name = 'C_D10N';
RunNeed.(char(name)) = {'openWnd','wrt_C_D10N'};
RunAxes.(char(name)) = struct('YLabel','',...
    'Title','Transfer coefficient for momentum at 10 m(neutral)');

% Transfer coefficient for heat at 10 m (neutral)
name = 'C_E10N';
RunNeed.(char(name)) = {'openWnd','wrt_C_E10N'};
RunAxes.(char(name)) = struct('YLabel','',...
    'Title','Transfer coefficient for heat at 10 m(neutral)');

% Transfer coefficient for humidity at 10 m (neutral)
name = 'C_H10N';
RunNeed.(char(name)) = {'openWnd','wrt_C_H10N'};
RunAxes.(char(name)) = struct('YLabel','',...
    'Title','Transfer coefficient for humidity at 10 m(neutral)');

% Transfer coefficient for momentum
name = 'C_D';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_C_D',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel','',...
    'Title','Transfer coefficient for momentum');

% Transfer coefficient for heat
name = 'C_E';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_C_E',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel','',...
    'Title','Transfer coefficient for heat');

% Transfer coefficient for humidity
name = 'C_H';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_C_H',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel','',...
    'Title','Transfer coefficient for humidity');

% Transfer coefficient for momentum at 10 m
name = 'C_D10';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_C_D10',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel','',...
    'Title','Transfer coefficient for momentum at 10 m');

% Transfer coefficient for heat at 10 m
name = 'C_E10';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_C_E10',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel','',...
    'Title','Transfer coefficient for heat at 10 m');

% Transfer coefficient for humidity at 10 m
name = 'C_H10';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_C_H10',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel','',...
    'Title','Transfer coefficient for humidity at 10 m');

% Wind speed at 10 m
name = 'u10';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_u10',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel',' (m s^{-1})',...
    'Title','Wind speed at 10 m');

% Wind speed at 10 m (neutral)
name = 'u10N';
RunNeed.(char(name)) = {'openWnd','wrt_u10N'};
RunAxes.(char(name)) = struct('YLabel',' (m s^{-1})',...
    'Title','Wind speed at 10 m (neutral)');

% Relative humidity at 10 m
name = 'rh10';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_rh10',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel',' (%%)',...
    'Title','Relative humidity at 10 m');

% Air temperature at 10 m
name = 't10';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_t10',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel',' (^o{C})',...
    'Title','Air temperature at 10 m');

% Long wave heat flux
name = 'Qlnet';
RunNeed.(char(name)) = {'openRH','openAirT','openWtr','openSW',...
    'openLWnet','wrt_Qlnet'};
RunAxes.(char(name)) = struct('YLabel',' (W m^{-2})',...
    'Title','Net long wave heat radiation');

% Incoming long wave radiation
name = 'Qlin';
RunNeed.(char(name)) = {'openRH','openAirT','openWtr','openSW',...
    'open_LWnet','wrt_Qlin'};
RunAxes.(char(name)) = struct('YLabel',' (W m^{-2})',...
    'Title','Incoming long wave radiation');

% Outgoing long wave heat flux
name = 'Qlout';
RunNeed.(char(name)) = {'openWtr','wrt_Qlout'};
RunAxes.(char(name)) = struct('YLabel',' (W m^{-2})',...
    'Title','Outgoing long wave radiation');

% Air shear velocity
name = 'uSt_a';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_uSt_a',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel',' (m s^{-1})',...
    'Title','Air shear velocity');

% Air shear velocity
name = 'uSt_aN';
RunNeed.(char(name)) = {'openWnd','wrt_uSt_aN'};
RunAxes.(char(name)) = struct('YLabel',' (m s^{-1})',...
    'Title','Air shear velocity (neutral)');

% Evaporation
name = 'Evap';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_Evap',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel',' (mm day^{-1})',...
    'Title','Evaporation');

% Reflected short wave radiation
name = 'Qsr';
RunNeed.(char(name)) = {'openSW','wrt_Qsr'};
RunAxes.(char(name)) = struct('YLabel',' (W m^{-2})',...
    'Title','Reflected short wave radiation');

% Monin-Obukhov length scale
name = 'obu';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','wrt_obu',...
    'senslatYes'};
RunAxes.(char(name)) = struct('YLabel',' ',...
    'Title','Atmospheric stability');

% total heat flux
name = 'Qtot';
RunNeed.(char(name)) = {'openWtr','openWnd','openRH','openAirT','openSW',...
    'openLWnet','wrt_Qtot','QtotYes'};
RunAxes.(char(name)) = struct('YLabel',' (W m^{-2})',...
    'Title','Total surface heat flux');

% Fill Construct
fNms = fieldnames(outputConstruct);
for k = 1:length(fNms)
    strN = RunNeed.(char(fNms{k}));
    for j = 1:length(strN)
        outputConstruct.(char(fNms{k})).(char(strN{j})) = true;
    end
    axN = fieldnames(RunAxes.(char(fNms{k})));
    for j = 1:length(axN)
        plotTable.(char(fNms{k})).(char(axN{j})) = ...
            RunAxes.(char(fNms{k})).(char(axN{j}));
    end
end

fNms = fieldnames(truthTable);
for k = 1:length(fNms)
    for j = 1:length(outputNames)
        if outputConstruct.(char(outputNames{j})).(char(fNms{k}))
            truthTable.(char(fNms{k})) = true;
        end
    end
end

fNms = fieldnames(writeTable);
for k = 1:length(fNms)
    for j = 1:length(outputNames)
        if outputConstruct.(char(outputNames{j})).(['wrt_' char(fNms{k})])
            writeTable.(char(fNms{k})) = true;
        end
    end
end
end

