%% Settings
%Simulation Settings
N =2; %Number of orbits
timesteps = 267; %Time steps per orbit (267 has the best results using Optimize_StepSize2VerifySTK_2Body.m)
Settings_General; %Loads paths
Settings_Radar; %Gives us wavelength and valid array spacing
Settings_Constellation; %Gives us const_struct
Settings_Targets; %Gives geocircles for targets

%% Propagate Orbits
%Matlab Simulation
[sats_matlab] = Sim_Constellation_matlab_2Body(sats,init_utcvec,N,timesteps);
%STK Simulation
[sats_stk] = Sim_Constellation_STK_2Body(sats,init_utcvec,N,timesteps);

%% Plot the comparison of matlab simulation and stk simulation
num_states = N*timesteps;
driving_sat_id = 1;%find(strcmp('TX',const_vec)==1,1,'first');
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
title(sprintf('Single Satellite Comparison\n Matlab vs. STK Simulation (Two-Body)'))

fprintf('Maximum Offset: %.f meters\n',max(magoff)*1000)