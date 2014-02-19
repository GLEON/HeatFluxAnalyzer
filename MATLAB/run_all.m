close all
clear all
TopFolder = 'C:\CEH Lancaster\Chapters\Bulk parameterization of the surface fluxes\HeatFlux Analyzer\MATLAB_v1';
cd(TopFolder);
lakeName = dir(TopFolder);
lakeName = {lakeName.name}; 
lakeName(strncmp(lakeName,'.',1)) = [];
lakeName(strncmp(lakeName,'Data',4)) = [];
lakeName = lakeName(1:38);
lakeName(27:28) = [];
lakeName([9,20,23]) = [];
lakeName(21) = [];

for i = 1:length(lakeName);
    Run_LHFA(lakeName{i},lakeName{i},true);
end