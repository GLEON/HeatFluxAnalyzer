function ...
    plotLA_results(writeTable,plotTable,dates,LakeName,Year)

fNms = fieldnames(writeTable);

figure_defaults = plotTable.FigD;
warning off all;
close all;
for k = 1:length(fNms)
    if ~islogical(writeTable.(char(fNms{k})))
        figure1 = figure;
        figFields = fieldnames(figure_defaults);
        for i = 1:length(figFields)
            set(figure1,figFields{i},figure_defaults.(char(figFields{i})))
        end
        axFields = fieldnames(plotTable.(char(fNms{k})));
        for i = 1:length(axFields)
            tempT = get(gca,axFields{i});
            if isequal(tempT,[0 1])
                set(gca,axFields{i},plotTable.(char(fNms{k})).(char(axFields{i})))
            elseif ishandle(tempT)
                set(tempT,'String',...
                    plotTable.(char(fNms{k})).(char(axFields{i})),...
                    'FontName',plotTable.(char(fNms{k})).FontName,...
                	'FontSize',plotTable.(char(fNms{k})).FontSize)
            else
                set(gca,axFields{i},plotTable.(char(fNms{k})).(char(axFields{i})))
            end
            hold on;
        end
        plot(dates,writeTable.(char(fNms{k})),'k','Parent',gca);
        datetick('x');
        set(gca,'XLim',[min(dates) max(dates)]);
        print(plotTable.PrintD.format,plotTable.PrintD.res,...
            [Year '/' LakeName '_' char(fNms{k})])
        if plotTable.PrintD.toClose
            close
        end
    end
end
warning on all;

