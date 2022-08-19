%% Settings
%Simulation Settings
N =2; %Number of orbits
timesteps = 206; %Time steps per orbit (206 has the best results using Optimize_StepSize2VerifySTK_2Body.m)

Settings_General; %Loads paths
Settings_Radar; %Gives us wavelength and valid array spacing
Settings_Constellation; %Gives us const_struct
Settings_Targets; %Gives geocircles for targets
Settings_STK; %Gives solar conditions

%% Propagate Orbits
%Matlab Simulation
[sats_matlab] = Sim_Constellation_matlab_J2Drag(sats,init_utcvec,N,timesteps);
%STK Simulation
[sats_stk] = Sim_Constellation_STK_J2Drag(sats,init_utcvec,N,timesteps);

%% Plot the comparison of matlab simulation and stk simulation
num_states = N*timesteps;
driving_sat_id = find(strcmp('TX',const_vec)==1,1,'first');
xoff = zeros(1,num_states); yoff = zeros(1,num_states); zoff = zeros(1,num_states); 
magoff = zeros(1,num_states);
tvec = zeros(1,num_states);
for stateid = 1:num_states
    xoff(stateid) = sats_matlab(driving_sat_id).states(stateid).R(1)-sats_stk(driving_sat_id).states(stateid).R(1);
    yoff(stateid) = sats_matlab(driving_sat_id).states(stateid).R(2)-sats_stk(driving_sat_id).states(stateid).R(2);
    zoff(stateid) = sats_matlab(driving_sat_id).states(stateid).R(3)-sats_stk(driving_sat_id).states(stateid).R(3);
    magoff(stateid) = norm([xoff(stateid) yoff(stateid) zoff(stateid)]);
    tvec(stateid) = sats_matlab(driving_sat_id).states(stateid).t;
end
figure(101)
plot([0:length(tvec)-1]./timesteps,xoff*1000); hold on
plot([0:length(tvec)-1]./timesteps,yoff*1000);
plot([0:length(tvec)-1]./timesteps,zoff*1000);
plot([0:length(tvec)-1]./timesteps,magoff*1000);
hold off
grid on
xlabel('Orbits (~)')
ylabel('Difference (Matlab-STK) [m]')
ylim([-100 100])
legend({'X','Y','Z','Total'})
title(sprintf('Single Satellite Comparison\n Matlab vs. STK Simulation (J2 + Drag)'))

fprintf('Maximum Offset: %.f meters @ %.0f timesteps\n',max(magoff)*1000,timesteps)