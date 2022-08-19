% clearvars


%% Define Monte Carlo Settings
%Define Monte Carlo parameters
if 1
    nruns = 5; %Random runs for every point
    pos_off_magvec = [0 0.5]/1000; %km
    vel_off_magvec = [0 0.05]/1000; %km/s
    npos = 20; %Number of position points
    nvel = 20; %Number of velocity points
    dirpath = 'MonteCarlo_Aug_18_micro';
elseif 1
    nruns = 5; %Random runs for every point
    pos_off_magvec = [0 2]/1000; %km
    vel_off_magvec = [0 0.5]/1000; %km/s
    npos = 20; %Number of position points
    nvel = 20; %Number of velocity points
    dirpath = 'MonteCarlo_Aug_18_small';
elseif 0
    nruns = 5; %Random runs for every point
    pos_off_magvec = [0 5]/1000; %km
    vel_off_magvec = [0 2]/1000; %km/s
    npos = 15; %Number of position points
    nvel = 15; %Number of velocity points
    dirpath = 'MonteCarlo_data_p15_v15_n5_Aug_12';
end

%Define things for messages
timebetweenupdates = 45; %seconds
timeperorbsim = 4.43; %sec/orbit Simulation run time estimate

%Define Orbit Simulation Settings
Ndays = 1; %365.25
Norbits = 2; %Orbits per run
timesteps = 206; %Time steps per orbit

%Define Save information

fnfrmt = @(p,v,n)sprintf('data_p%.0f_v%.0f_n%.0f.mat',p,v,n);

savedatapath = fullfile('G:\Dissertation_Data',dirpath);
% savepath = fullfile(cd,'Results/SupportingData');
% savedatapath = fullfile(savepath,dirpath);

%% Load Settings
Settings_General; %Loads paths
Settings_Radar; %Gives us wavelength and valid array spacing
Settings_Constellation; %Gives us const_struct
Settings_Targets; %Gives geocircles for targets
Settings_STK

%% Precalculations
num_sats = length(sats);
%Axis vectors
pos_off_vec = linspace(min(pos_off_magvec),max(pos_off_magvec),npos);
vel_off_vec = linspace(min(vel_off_magvec),max(vel_off_magvec),nvel);


%% Monte Carlo Simulation
run_montecarlo = false; rerun=false;
totalruns = npos*nvel*nruns;
if ~exist(savedatapath,'dir') || run_montecarlo
    %Print the prep info
    totalorbitssimmed = totalruns*Norbits;
    timeestmin = timeperorbsim*totalorbitssimmed/60; %min
    if timeestmin < 60
        timeest = timeestmin; timeeststr = 'min';
    else
        timeest = timeestmin/60; timeeststr = 'hrs';
    end
    cmdsize = matlab.desktop.commandwindow.size;
    cmdline = repmat('=',1,cmdsize(1));
    fprintf('\n%s\n',cmdline);
    fprintf('This script will simulate %.0f constellation orbits with an estimated time of %.2f %s\n',...
        totalorbitssimmed,timeest,timeeststr);
    fprintf('%s\n\n',cmdline);
    
    %Confirm file overwrite
    if 1 
        if exist(savedatapath,'dir')
            if rerun
              fprintf('Directory %s will be overwritten.\n\tWould you like to continue? (y/n)',dirpath);
              runconfirm = input('','s');
              if ~strcmp(runconfirm,'y')
                  fprintf('Simulation Cancelled\n')
                  return
              end
            else
              fprintf('Running remaining jobs in %s\n',dirpath)
            end
        else
            %Make new directory
            mkdir(savedatapath)
            fprintf('Directory %s was created.\n',dirpath);
        end
    end
    if 1
        %Confirm run
        fprintf('Are you ready for this simulation run? (y/n)')
        runconfirm = input('','s');
        if strcmp(runconfirm,'n')
            fprintf('Simulation Cancelled\n')
            return
        else
            fprintf('Let''s a GO!\n')
        end
    end
    
    %Ensure parallel pool processing is off
    if ~isempty(gcp('nocreate'))
        %Delete the pool
        delete(gcp('nocreate'));
    end
    %Grab the cluster
    clust = parcluster();
    %Delete any lingering jobs
    delete(clust.Jobs);
    jobcheck = false; %true indicates valid job has been created; false means one will be created
    maxtasks = 40; %Maximum tasks 
    waittic = 5; %seconds
    timebetweenupdates = 60; %seconds
    
    
    %Iterate over position error vector
    sats_ideal = sats;
    total_n = 0;
    job_n = 0;
    totaljobs = ceil(totalruns/maxtasks);
    tic;
    %Load tasks to job
    for pid = 1:npos
        pos_off_rad = pos_off_vec(pid);
        %Iterate over velocity error vector
        for vid = 1:nvel
            vel_off_rad = vel_off_vec(vid);
            for n = 1:nruns
                total_n = total_n+1;
                sats_now = sats_ideal;
                if 1
                    if ~jobcheck
                        job_n = job_n+1;
                        %Send message
                        elapstr = parsetoc(toc);
                        fprintf('Creating Job %.0f of %.0f\n\t%s\n',job_n,totaljobs,elapstr)
                        %Create a new job
                        job = createJob(clust,'AttachedFiles','Subfunctions');
                        jobcheck = true;
                    end
                    %Check if task was previously completed
                    if exist(fullfile(savedatapath,fnfrmt(pid,vid,n)),'file')
                        %Skip this task
                        fprintf('File %s already exists, skipping...\n',fnfrmt(pid,vid,n))
                    else
                        %Send this task to the job
                        createTask(job,@MonteCarlo_func,1,...
                            {sats_now,pos_off_rad,vel_off_rad,init_utcvec,...
                            Norbits,timesteps,savedatapath,fnfrmt,pid,vid,n});
                    end

                    %Check for number of tasks and run at maxtasks
                    if length(job.Tasks) >= maxtasks
                        run_job_script
                        jobcheck = false;
                    end
                else
                    MonteCarlo_func(sats_now,pos_off_rad,vel_off_rad,init_utcvec,...
                        Norbits,timesteps,savedatapath,fnfrmt,pid,vid,n);
                end
            end
        end
    end
    %Ensure job is empty or run final tasks
    if jobcheck
        %Run the job
        run_job_script
        jobcheck = false;
    end
    
    %Send finishing message
    lasttoc = toc;
    fprintf('\n%s\n',cmdline);
    fprintf('This script has completed %0.f orbits in %.2f min\n\t Approximately %.2f sec/orbit\n',...
        totalorbitssimmed,lasttoc/60,lasttoc/totalorbitssimmed);
    
else
    %% Load propagation data
%     load(savedatapath)
end
%% Post Calculations
%Ensure parallel processing is on
if isempty(gcp('nocreate'))
    %Turn on pool
    parpool_now = parpool;
else
    %Grab pool
    parpool_now = gcp('nocreate');
end
%Initialize output variables
fails_out = zeros(npos,nvel);
fails_avg_time = zeros(npos,nvel);
fails_max_time = zeros(npos,nvel);
data_avg_nofail_out = zeros(npos,nvel);
data_max_nofail_out = zeros(npos,nvel);
truedata_avg_out = zeros(npos,nvel);
truedata_max_out = zeros(npos,nvel);
dataUML_avg_out = zeros(npos,nvel);
dataUML_max_out = zeros(npos,nvel);
trueUML_avg_out = zeros(npos,nvel);
trueUML_max_out = zeros(npos,nvel);

UML_cells = cell(npos,nvel);
data_cells = cell(npos,nvel);
failcheck_cells = cell(npos,nvel,nruns);
total_n = 0;
tic; 
lasttoc = toc;
for pid = 1:npos
    for vid = 1:nvel        
        %Initialize main outputs
        fails_n = 0; data_collect_nofail_n = zeros(1,nruns); dataUML_n = zeros(1,nruns); failUML_n = zeros(1,nruns);
        trueUML_n = zeros(1,nruns); truedata_n = zeros(1,nruns);
        %Send message            
        if toc-lasttoc > timebetweenupdates || (pid==1 && vid==1)
            %Print update
            fprintf('Constellation Calcs %.0f/%.0f:\n\tpos %.0f/%.0f, vel %.0f/%.0f\n\tTime Elapsed: %.2f min.\n',...
                total_n,totalruns,pid,npos,vid,nvel,toc/60);
            %Update last toc
            lasttoc = toc;
        end
        total_n = total_n+nruns;
        %Do output calculations
        parfor nid = 1:nruns
            %Load propagation data
            datafn = fullfile(savedatapath,feval(fnfrmt,pid,vid,nid));
            sats_now = dload(datafn);
%             sats_now = datanow.sats_now;
%             sats_now = sats_cells{pid,vid,nid};
            %Do constellation calculations
            [const_calc,failcheck_mat] = constellation_calcs(sats_now,driving_sat_id,avg_array_rangewl*wavelength,fail_min_off,fail_prob);
            
            %Deterime UML based on failure
            firstfailind = find(any(failcheck_mat'),1,'first');
            if isempty(firstfailind)
                firstfailind = length(const_calc.data_total);
            end
            failUML_n(nid) = sats_now(driving_sat_id).states(firstfailind).t;
            
            %Determine UML based on data
            lastdatind = find(diff(const_calc.data_total)>0,1,'last');
            if isempty(lastdatind)
                lastdatind = 1;
            end
            dataUML_n(nid) = sats_now(driving_sat_id).states(lastdatind).t;
            
            %Determine true UML based on failure and data
            trueendind = min([firstfailind lastdatind]);
            trueUML_n(nid) = sats_now(driving_sat_id).states(trueendind).t;
            
            %Send details to output variable
            data_collect_nofail_n(nid) = const_calc.data_total(lastdatind);
            truedata_n(nid) = const_calc.data_total(trueendind);
            
            %Send out save data
            failcheck_cells{pid,vid,nid} = failcheck_mat;
%             const_calc_cells{pid,vid,nid} = const_calc;
            if any(any(failcheck_cells{pid,vid,nid}))
                fails_n = fails_n+1;
            end
            
        end
        
        fails_out(pid,vid) = fails_n;
        data_avg_nofail_out(pid,vid) = mean(data_collect_nofail_n);
        data_max_nofail_out(pid,vid) = max(data_collect_nofail_n);
        truedata_avg_out(pid,vid) = mean(truedata_n);
        truedata_max_out(pid,vid) = max(truedata_n);
        dataUML_avg_out(pid,vid) = mean(dataUML_n);
        dataUML_max_out(pid,vid) = max(dataUML_n);
        trueUML_avg_out(pid,vid) = mean(trueUML_n);
        trueUML_max_out(pid,vid) = max(trueUML_n);
    end
end

%% Plot the results
fignum = 500;

%Data UML Average vs. each accuracy point
% fignum = fignum+1;
% figure(fignum); clf
% contourf(pos_off_vec*1000,vel_off_vec*1000,dataUML_avg_out')
% grid on
% xlabel('Position Accuracy (m)')
% ylabel('Velocity Accuracy (m)')
% c = colorbar;
% c.Label.String = 'Average UML (sec)';
% title(sprintf('Data Useful Mission Lifetime Average'))
% FontWidthandPos

%Data UML Max vs. each accuracy point
% fignum = fignum+1;
% figure(fignum); clf
% contourf(pos_off_vec*1000,vel_off_vec*1000,dataUML_max_out')
% xlabel('Position Accuracy (m)')
% ylabel('Velocity Accuracy (m)')
% c = colorbar;
% c.Label.String = 'Maximum UML (sec)';
% title(sprintf('Data Useful Mission Lifetime Maximum'))
% FontWidthandPos

%Failures vs. each accuracy point
fignum = fignum+1;
figure(fignum); clf
contourf(pos_off_vec*1000,vel_off_vec*1000,fails_out')
xlabel('Position Accuracy (m)')
ylabel('Velocity Accuracy (m)')
c = colorbar;
c.Label.String = '# of Failures';
title('Failures')
grid on 
FontWidthandPos

%True UML Average vs. each accuracy point
fignum = fignum+1;
figure(fignum); clf
contourf(pos_off_vec*1000,vel_off_vec*1000,trueUML_avg_out')
xlabel('Position Accuracy (m)')
ylabel('Velocity Accuracy (m)')
c = colorbar;
c.Label.String = 'Average UML (sec)';
title(sprintf('Useful Mission Lifetime Average'))
grid on 
FontWidthandPos

%True UML Max vs. each accuracy point
% fignum = fignum+1;
% figure(fignum); clf
% contourf(pos_off_vec*1000,vel_off_vec*1000,trueUML_max_out')
% xlabel('Position Accuracy (m)')
% ylabel('Velocity Accuracy (m)')
% c = colorbar;
% c.Label.String = 'Maximum UML (sec)';
% title(sprintf('Useful Mission Lifetime Maximum'))

%Average data (nofail) vs. each accuracy point
% fignum = fignum+1;
% figure(fignum); clf
% contourf(pos_off_vec*1000,vel_off_vec*1000,data_avg_nofail_out')
% xlabel('Position Accuracy (m)')
% ylabel('Velocity Accuracy (m)')
% c = colorbar;
% c.Label.String = 'Average Data Collected (line-km)';
% title(sprintf('Average Data Collected\nIgnoring Failures'))

%Max data (nofail) vs. each accuracy point
% fignum = fignum+1;
% figure(fignum); clf
% contourf(pos_off_vec*1000,vel_off_vec*1000,data_max_nofail_out')
% xlabel('Position Accuracy (m)')
% ylabel('Velocity Accuracy (m)')
% c = colorbar;
% c.Label.String = 'Maximum Data Collected (line-km)';
% title(sprintf('Maximum Data Collected\nIgnoring Failures'))

%Average data (true) vs. each accuracy point
fignum = fignum+1;
figure(fignum); clf
contourf(pos_off_vec*1000,vel_off_vec*1000,truedata_avg_out')
xlabel('Position Accuracy (m)')
ylabel('Velocity Accuracy (m)')
c = colorbar;
c.Label.String = 'Average Data Collected (line-km)';
title(sprintf('Average Data Collected'))
grid on
FontWidthandPos

%Max data (true) vs. each accuracy point
% fignum = fignum+1;
% figure(fignum); clf
% contourf(pos_off_vec*1000,vel_off_vec*1000,truedata_max_out')
% xlabel('Position Accuracy (m)')
% ylabel('Velocity Accuracy (m)')
% c = colorbar;
% c.Label.String = 'Maximum Data Collected (line-km)';
% title(sprintf('Maximum Data Collected\nConsidering Failures'))