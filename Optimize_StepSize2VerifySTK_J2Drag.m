tic
timestepsvec = [160:300];%;
N = 2;

%Initialize Settings
Settings_General; %Loads paths
Settings_Radar; %Gives us wavelength and valid array spacing
Settings_Constellation; %Gives us const_struct
Settings_Targets; %Gives geocircles for targets

magout = zeros(size(timestepsvec));
%Iterate
for tid = 1:length(timestepsvec)
    timesteps = timestepsvec(tid);

    %% Propagate Orbits
    %Matlab Simulation
    [sats_matlab] = Sim_Constellation_matlab_J2Drag(sats,init_utcvec,N,timesteps);
    %STK Simulation
    [sats_stk] = Sim_Constellation_STK_J2Drag(sats,init_utcvec,N,timesteps);

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
    fprintf('Maximum Offset: %.f meters @ %.0f timesteps\n',max(magoff)*1000,timesteps)
    magout(tid) = max(magoff);
end
figure(1001)
plot(timestepsvec,magout*1000)
xlabel('Time Steps Per Orbit')
ylabel('Maximum Offset (Matlab-STK)')
toc