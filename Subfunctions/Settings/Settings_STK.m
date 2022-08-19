%% Settings
%STK Propagation Settings
statetype = 'Cartesian'; %'Classical'
propagatortype = 'HPOP'; %'TwoBody'
coordsyststr = 'J2000'; %'J2000' 'ICRF'
% daystoend = 30; %a month

%File details
scen_name = 'Dissertation_Sim_v1';
savedir = fullfile(cd,'STK_Scenarios');
savepath_stk = fullfile(savedir,[scen_name '.sc']);

%Model details
modelbasedir = 'C:\Program Files\AGI\STK 12\STKData\VO\Models\Space';
RXmodelfn = 'cubesat_3U.dae'; 
TXmodelfn = 'cubesat_6U.dae';

%Drag solar settings
F10abs = 78;
F10avg = 78.71;
Apvec = [3 3 3 2 3 3 4];
%% Automated Section

