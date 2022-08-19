clc
clearvars

addpath(genpath('Subfunctions'))
%% Settings
%Simulation Settings
% stepsforward = 150;
Norbits =120; %Number of orbits
timesteps = 206; %Time steps per orbit
% N = stepsforward/timesteps;

Settings_General; %Loads paths
Settings_Radar; %Gives us wavelength and valid array spacing
Settings_Constellation; %Gives us const_struct
% Settings_Targets; %Gives geocircles for targets
% Settings_STK


%Set allowables for failure
fail_dist_m = 0.25; %m Distance that satellites must be separated by at all times
fail_prob = 1e8; %Always fail

%Define things for messages
timebetweenupdates = 45; %seconds
timeperorbsim = 4.43; %sec/orbit Simulation run time estimate

%Define things for cluster
maxtasks = 40; %Maximum tasks 
waittic = 1; %seconds

%% Define boundaries

%Define offsets of interest and save information
if 1
  %This is the big sweep result for the dissertation (note the best data
  %and best UML result) %Do the nu of 5 steps to finish this
  Norbits =30; %Number of orbits
    wsteps = 5;
    woff_min = .000002*0;
    woff_max = .0000002*1; 
    nusteps = 5;
    nuoff_min = 0.000006*0;
    nuoff_max = 1e-7*1;
    isteps = 10;
    ioff_min = 0.5e-7;
    ioff_max = 1.5e-6;
    RAsteps = 5;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.0000005*1;
    dirpath = 'Optimize_8_15_Big_Run_Take_3';
elseif 1
%This is the big sweep result for the dissertation (note the best data
  %and best UML result) %Do the nu of 5 steps to finish this
  Norbits =30; %Number of orbits
    wsteps = 1;
    woff_min = .000002*0;
    woff_max = .0000002*0; 
    nusteps = 10;
    nuoff_min = 0.000006*0;
    nuoff_max = 1e-7*1;
    isteps = 10;
    ioff_min = 10e-7;
    ioff_max = 30e-6;
    RAsteps = 10;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.0000005*1;
    dirpath = 'Optimize_8_15_inc_vs_RA';
elseif 1
  %This is the peek result only looking at optimizing UML baed on the big
  %sweep result
    wsteps = 5;
    woff_min = .000002*0;
    woff_max = .00000075*1; 
    nusteps = 5;
    nuoff_min = 0.000006*0;
    nuoff_max = 0.0000005*1;
    isteps = 9;
    ioff_min = 1e-6;
    ioff_max = 5e-6;
    RAsteps = 5;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.0000015*1;
    dirpath = 'Optimize_8_12_UML_Sweep_high_inc';
elseif 0
  %This is the peek result only looking at optimizing UML baed on the big
  %sweep result
    wsteps = 1;
    woff_min = .000002*0;
    woff_max = .00000075*0; 
    nusteps = 1;
    nuoff_min = 0.000006*0;
    nuoff_max = 0.0000005*0;
    isteps = 40;
    ioff_min = 3.2e-6;
    ioff_max = 5e-6;
    RAsteps = 1;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.0000015*0;
    dirpath = 'Optimize_8_12_UML_Peek_inc_only_wide';
elseif 0
  %This is the peek result only looking at optimizing UML baed on the big
  %sweep result
    wsteps = 1;
    woff_min = .000002*0;
    woff_max = .00000075*0; 
    nusteps = 1;
    nuoff_min = 0.000006*0;
    nuoff_max = 0.0000005*0;
    isteps = 50;
    ioff_min = 1.2e-6;
    ioff_max = 3.2e-6;
    RAsteps = 1;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.0000015*0;
    dirpath = 'Optimize_8_12_UML_Peek_inc_only';
elseif 1
  %This is the big sweep result for the dissertation (note the best data
  %and best UML result) %Do the nu of 5 steps to finish this
    wsteps = 5;
    woff_min = .000002*0;
    woff_max = .00000075*1; 
    nusteps = 5;
    nuoff_min = 0.000006*0;
    nuoff_max = 0.0000005*1;
    isteps = 10;
    ioff_min = 0.5e-7;
    ioff_max = 1.5e-6;
    RAsteps = 5;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.0000015*1;
    dirpath = 'Optimize_8_10_Big_Run';
elseif 1
    wsteps = 5;
    woff_min = .000002*0;
    woff_max = .000002*1; 
    nusteps = 1;
    nuoff_min = 0.000006*0;
    nuoff_max = 0.000006*0;
    isteps = 10;
    ioff_min = 2e-7;
    ioff_max = 2e-6;
    RAsteps = 1;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.000001*0;
    dirpath = 'Optimize_8_8_w';
elseif 0
    wsteps = 1;
    woff_min = .000002*0;
    woff_max = .000002*0; 
    nusteps = 5;
    nuoff_min = 0.000006*0;
    nuoff_max = 0.000006*1;
    isteps = 10;
    ioff_min = 2e-7;
    ioff_max = 2e-6;
    RAsteps = 1;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.000001*0;
    dirpath = 'Optimize_8_8_nu';
elseif 1
    wsteps = 1;
    woff_min = .000002*0;
    woff_max = .000002*0; 
    nusteps = 1;
    nuoff_min = 0.000006*0;
    nuoff_max = 0.000006*0;
    isteps = 10;
    ioff_min = 2e-7;
    ioff_max = 2e-6;
    RAsteps = 5;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.000001*1;
    dirpath = 'Optimize_8_8_RA';
elseif 0 %Good Inclination profile (use in report)
    wsteps = 1;
    woff_min = .000002*0;
    woff_max = .000002*0; 
    nusteps = 1;
    nuoff_min = 0.000006*0;
    nuoff_max = 0.000006*0;
    isteps = 21;
    ioff_min = 0.000001*0;
    ioff_max = 0.0000002*20;
    RAsteps = 1;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.000001*0;
    dirpath = 'Optimize_8_8';
elseif 0
     wsteps = 2;
    woff_min = .000002*0;
    woff_max = .0000002*2; 
    nusteps = 2;
    nuoff_min = 0.000006*0;
    nuoff_max = 0.0000006*1;
    isteps = 40;
    ioff_min = 0.000001*0;
    ioff_max = 0.0000001*6;
    RAsteps = 2;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.0000001*2;
    dirpath = 'Optimize_Incline_7_18';
elseif 1
    wsteps = 5;
    woff_min = .000002*0;
    woff_max = .000002*4; 
    nusteps = 3;
    nuoff_min = 0.000006*0;
    nuoff_max = 0.000006*2;
    isteps = 21;
    ioff_min = 0.000001*0;
    ioff_max = 0.000001*20;
    RAsteps = 5;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.000001*4;
    dirpath = 'Optimize_7_15';
elseif 1
    wsteps = 10;
    woff_min = .000002*4;
    woff_max = .000002*8; 
    nusteps = 4;
    nuoff_min = 0.000006*2;
    nuoff_max = 0.000006*4;
    isteps = 10;
    ioff_min = 0.000001*12;
    ioff_max = 0.000001*24;
    RAsteps = 4;
    RAoff_min = 0.000001*4;
    RAoff_max = 0.000001*8;
    dirpath = 'Optimize_7_13_2';
elseif 1
    wsteps = 10;
    woff_min = .000002*2;
    woff_max = .000002*4; 
    nusteps = 4;
    nuoff_min = 0.000006*1;
    nuoff_max = 0.000006*2;
    isteps = 10;
    ioff_min = 0.000001*6;
    ioff_max = 0.000001*12;
    RAsteps = 4;
    RAoff_min = 0.000001*2;
    RAoff_max = 0.000001*4;
    dirpath = 'Optimize_7_13';
elseif 1
    wsteps = 10;
    woff_min = .000002*0;
    woff_max = .000002*2; 
    nusteps = 4;
    nuoff_min = 0.000006*0;
    nuoff_max = 0.000006*1;
    isteps = 10;
    ioff_min = 0.000001*0;
    ioff_max = 0.000001*6;
    RAsteps = 4;
    RAoff_min = 0.000001*0;
    RAoff_max = 0.000001*2;
    dirpath = 'Optimize_7_8';
end
    
% savedatapath = fullfile(cd,'Results/SupportingData',dirpath);
savedatapath = fullfile('G:\Dissertation\data',dirpath);
fnfrmt = @(wid,nuid,iid,raid)sprintf('data_w%.0f_nu%.0f_i%.0f_RA%.0f.mat',wid,nuid,iid,raid);

%% Optimize Loop
%Bound independent variables
woffs = linspace(woff_min,woff_max,wsteps);
nuoffs = linspace(nuoff_min,nuoff_max,nusteps);
ioffs = linspace(ioff_min,ioff_max,isteps);
RAoffs = linspace(RAoff_min,RAoff_max,RAsteps);

run_optimize = false; rerun = false;
totalruns = wsteps*nusteps*isteps*RAsteps;
totalorbitssimmed = totalruns*Norbits;

if ~exist(savedatapath,'dir') || run_optimize
    %Print the prep info
    
    timeestmin = timeperorbsim*totalorbitssimmed/60; %min
    if timeestmin < 60
        timeest = timeestmin; timeeststr = 'min';
    else
        timeest = timeestmin/60; timeeststr = 'hrs';
    end
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
    
    %Necessary Precalculations
    avg_array_range_m = avg_array_rangewl*wavelength;

    ntotal = wsteps*nusteps*isteps*RAsteps;
    
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
 
    %Iterate
    total_n = 0;
    job_n = 0;
    totaljobs = ceil(totalruns/maxtasks);
    total_n = 0; 
    tic;
    lasttoc = toc;
    timebetweenupdates = 30;
    est_dt_n = 3; %sec
    for wid = 1:wsteps
        wangs = woffs(wid)*((1:numel(const_vec))-1);
        for nuid = 1:nusteps
            nuangs = nuoffs(nuid)*((1:numel(const_vec))-1);
            for iid = 1:isteps
                iangs = ioffs(iid)*((1:numel(const_vec))-1);
                for RAid = 1:RAsteps
                    RAangs = RAoffs(RAid)*((1:numel(const_vec))-1);
                    
                    %Increase run count
                    total_n = total_n+1;
                    
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
                        if exist(fullfile(savedatapath,fnfrmt(wid,nuid,iid,RAid)),'file')
                            %Skip this task
                            fprintf('File %s already exists, skipping...\n',fnfrmt(wid,nuid,iid,RAid))
                        else
                            %Send this task to the job
                            createTask(job,@optimize_func,0,...
                                {const_vec,OE0,nuangs,iangs,RAangs,wangs,init_utcvec,...
                                Norbits,timesteps,avg_array_range_m,fail_dist_m,fail_prob,...
                                savedatapath,fnfrmt,wid,nuid,iid,RAid});
                        end

                        %Check for number of tasks and run at maxtasks
                        if length(job.Tasks) >= maxtasks
                            run_job_script
                            jobcheck = false;
                        end
                    else
                        optimize_func(const_vec,OE0,nuangs,iangs,RAangs,wangs,init_utcvec,...
                            Norbits,timesteps,avg_array_range_m,fail_dist_m,fail_prob,...
                            savedatapath,fnfrmt,wid,nuid,iid,RAid);
                    end
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
    tic
end
%% Load the data
%Initialize outputs
dataout = zeros(wsteps,nusteps,isteps,RAsteps);
failout = zeros(wsteps,nusteps,isteps,RAsteps);
minoffout = zeros(wsteps,nusteps,isteps,RAsteps);
mintout = zeros(wsteps,nusteps,isteps,RAsteps);
lid = 0;
for wid = 1:wsteps
    for nuid = 1:nusteps
        for iid = 1:isteps
            for RAid = 1:RAsteps
                %Print Message 
                lid = 1 + lid;
                fprintf('Loading %d of %d data files\n',lid,numel(dataout));
                %Load save data
                savefn = feval(fnfrmt,wid,nuid,iid,RAid);
                datanow = load(fullfile(savedatapath,savefn));
                %Find if failed
                failout(wid,nuid,iid,RAid) = any(any(datanow.failcheck_mat));
                %Find when it failed
                UMLid = find(any(datanow.failcheck_mat')==1,1,'first');
                if isempty(UMLid)
                    UMLid = size(datanow.failcheck_mat,1);
                end
                UMLout(wid,nuid,iid,RAid) = UMLid/timesteps; %In number of orbits
                %Find total data collected
                dataout(wid,nuid,iid,RAid) = datanow.out_struct.data_total(UMLid);
                %Find minimum distance and time
                minoffout(wid,nuid,iid,RAid) = datanow.minnow;
                mintout(wid,nuid,iid,RAid) = datanow.mintnow;
            end
        end
    end
end

lasttoc = toc;
fprintf('\n%s\n',cmdline);
fprintf('This script has loaded %0.f orbits in %.2f min\n\t Approximately %.2f sec/orbit\n',...
    totalorbitssimmed,lasttoc/60,lasttoc/totalorbitssimmed);
%% Plot the single line trend
% %Plot the wang
% figure(111);clf
% subplot(1,4,1)
% plot(woffs,reshape(dataout(:,1,1,1),1,wsteps));
% xlabel('W Off (rad)')
% ylabel('Data (line-km)')
% title('Data Out')
% subplot(1,4,2)
% plot(woffs,reshape(failout(:,1,1,1),1,wsteps));
% xlabel('W Off (rad)')
% ylabel('1=Fail, 0=Success')
% title('Fail Out')
% subplot(1,4,3)
% plot(woffs,reshape(minoffout(:,1,1,1),1,wsteps));
% xlabel('W Off (rad)')
% ylabel('Min. Offset (m)')
% title('Min. Offset Out')
% subplot(1,4,4)
% plot(woffs,reshape(mintout(:,1,1,1),1,wsteps));
% xlabel('W Off (rad)')
% ylabel('Min. Offset (sec)')
% title('Min. Time Out')
% 
% %Plot the nuang
% figure(222);clf
% subplot(1,4,1)
% plot(nuoffs,reshape(dataout(1,:,1,1),1,nusteps));
% xlabel('Nu Off (rad)')
% ylabel('Data (line-km)')
% title('Data Out')
% subplot(1,4,2)
% plot(nuoffs,reshape(failout(1,:,1,1),1,nusteps));
% xlabel('Nu Off (rad)')
% ylabel('1=Fail, 0=Success')
% title('Fail Out')
% subplot(1,4,3)
% plot(nuoffs,reshape(minoffout(1,:,1,1),1,nusteps));
% xlabel('Nu Off (rad)')
% ylabel('Min. Offset (m)')
% title('Min. Offset Out')
% subplot(1,4,4)
% plot(nuoffs,reshape(mintout(1,:,1,1),1,nusteps));
% xlabel('Nu Off (rad)')
% ylabel('Min. Offset (sec)')
% title('Min. Time Out')

%Plot the iang for each sub element
fignum = 333;
leg_pres = {'w','\nu','\Omega'};
plot_vecs = {woffs,nuoffs,RAoffs};
for pid = 1:length(leg_pres)
    fignum = fignum+1;
    evec = plot_vecs{pid};
    %Initiate Plotting
    figure(fignum);clf
    legs = {};
    
    if 1 && length(evec)>1
      %Plot as contour plots
      %Get Plot indices
      basesubs = cell(1,5);
      basesubs{1} = size(dataout);
      isubs = 1:size(dataout,3);
      switch pid
        case 1
          wsubs = 1:size(dataout,1);
          [imat,wmat] = meshgrid(isubs,wsubs);
          omat = ones(size(imat));
          basesubs = {size(dataout),wmat,omat,imat,omat};
          plotinds = sub2ind(basesubs{:});
        case 2
          nusubs = 1:size(dataout,2);
          [imat,numat] = meshgrid(isubs,nusubs);
          omat = ones(size(imat));
          basesubs = {size(dataout),omat,numat,imat,omat};
          plotinds = sub2ind(basesubs{:});
        case 3
          RAsubs = 1:size(dataout,4);
          [imat,RAmat] = meshgrid(isubs,RAsubs);
          omat = ones(size(imat));
          basesubs = {size(dataout),omat,omat,imat,RAmat};
          plotinds = sub2ind(basesubs{:});
      end
      %Do the plotting
      subplot(1,3,1)
      contourf(ioffs,evec,reshape(dataout(plotinds),length(evec),isteps),20);
      subplot(1,3,2)
      contourf(ioffs,evec,reshape(UMLout(plotinds),length(evec),isteps),20);
      subplot(1,3,3)
      contourf(ioffs,evec,reshape(minoffout(plotinds),length(evec),isteps),20);
%       %Do formatting
      subplot(1,3,1)
      hold off
      grid on
      xlabel('i Off (rad)')
      ylabel(sprintf('%s off (rad)',leg_pres{pid}))
      c = colorbar;
      c.Label.String = 'Data (line-km)';
      caxis([1e4 10e4])
      title('Data Out')
      subplot(1,3,2)
      hold off
      grid on
      xlabel('i Off (rad)')
      ylabel(sprintf('%s off (rad)',leg_pres{pid}))
      c = colorbar;
      c.Label.String = 'UML # of Orbits';
      caxis([2 20])
      title('UML')
      subplot(1,3,3)
      hold off
      grid on
      xlabel('i Off (rad)')
      ylabel(sprintf('%s off (rad)',leg_pres{pid}))
      c = colorbar;
      c.Label.String = 'Min. Offset (m)';
      caxis([0.2 0.3])
      title('Min. Offset Out')
      
      posvec = [165 508 1703 389];%[275 735 2278 522];
      FontWidthandPos(posvec)
    elseif 0
      %Plot as individual lines
      for eid = 1:length(evec)
          %Make legend entry
          legs{end+1} = sprintf('%s=%.2e rad.',leg_pres{pid},evec(eid));
          %Make cell array for grabbing plotting points
          basesubs = {size(dataout),ones(1,isteps),ones(1,isteps),1:size(dataout,3),ones(1,isteps)};
          if pid <3
              basesubs{pid+1} = basesubs{pid+1}*eid;
          else
              basesubs{pid+2} = basesubs{pid+2}*eid;
          end
          plotinds = sub2ind(basesubs{:});

          %Do the plotting
          subplot(1,4,1)

          plot(ioffs,reshape(dataout(plotinds),1,isteps));
          hold on
          subplot(1,4,2)
          plot(ioffs,reshape(UMLout(plotinds),1,isteps));
          hold on
          subplot(1,4,3)
          plot(ioffs,reshape(minoffout(plotinds),1,isteps));
          hold on
          subplot(1,4,4)
          plot(ioffs,reshape(mintout(plotinds),1,isteps));
          hold on
      end
      %Do formatting
      subplot(1,4,1)
      legend(legs)
      hold off
      grid on
      xlabel('i Off (rad)')
      ylabel('Data (line-km)')
      title('Data Out')
      subplot(1,4,2)
      legend(legs)
      hold off
      grid on
      xlabel('i Off (rad)')
      ylabel('UML (# of Orbits)')
      title('UML')
      subplot(1,4,3)
      legend(legs)
      hold off
      grid on
      xlabel('i Off (rad)')
      ylabel('Min. Offset (m)')
      title('Min. Offset Out')
      subplot(1,4,4)
      legend(legs)
      hold off
      grid on
      xlabel('i Off (rad)')
      ylabel('Min. Offset (sec)')
      title('Min. Time Out')
    end
end
% 
% %Plot the raang
% figure(444);clf
% subplot(1,4,1)
% plot(RAoffs,reshape(dataout(1,1,1,:),1,RAsteps));
% xlabel('RA Off (rad)')
% ylabel('Data (line-km)')
% title('Data Out')
% subplot(1,4,2)
% plot(RAoffs,reshape(failout(1,1,1,:),1,RAsteps));
% xlabel('RA Off (rad)')
% ylabel('1=Fail, 0=Success')
% title('Fail Out')
% subplot(1,4,3)
% plot(RAoffs,reshape(minoffout(1,1,1,:),1,RAsteps));
% xlabel('RA Off (rad)')
% ylabel('Min. Offset (m)')
% title('Min. Offset Out')
% subplot(1,4,4)
% plot(RAoffs,reshape(mintout(1,1,1,:),1,RAsteps));
% xlabel('RA Off (rad)')
% ylabel('Min. Offset (sec)')
% title('Min. Time Out')

%% Print best result
%Find maximum data collected without failure
maxdata = max(max(max(max(dataout))));
maxmindist = max(max(max(max(minoffout))));
maxUML = max(max(max(max(UMLout))));
%Find indices of best data point
bestdataind = find(dataout==maxdata);
[wsub,nusub,isub,rasub] = ind2sub(size(dataout),bestdataind);
% bestdataind = find(dataout==maxmindist);
% [wsub,nusub,isub,rasub] = ind2sub(size(dataout),bestdataind);
%Grab the bestdata values
wbestdata = woffs(wsub); nubestdata = nuoffs(nusub); ibestdata = ioffs(isub); RAbestdata = RAoffs(rasub);
databestdata = dataout(wsub,nusub,isub,rasub);
UMLbestdata = UMLout(wsub,nusub,isub,rasub);
minoffbestdata = minoffout(wsub,nusub,isub,rasub);
mintbestdata = mintout(wsub,nusub,isub,rasub);
%Print outputs
fprintf('Max Data Result: %.0e line-km\nUML: %.2f Orbits\nMin. Offset: %.2f m\nTime of Min.: %.0f sec\nwoff: %.3e rad.\nnuoff: %.3e rad.\nioff: %.3e rad.\nRAoff: %.3e rad.\n\n',...
    databestdata,UMLbestdata,minoffbestdata,mintbestdata,wbestdata,nubestdata,ibestdata,RAbestdata)
  
%Find inidices of best UML point
bestUMLind = find(UMLout==maxUML);
[wsub,nusub,isub,rasub] = ind2sub(size(dataout),bestUMLind(1));
%Grab the bestUML values
wbestUML = woffs(wsub); nubestUML = nuoffs(nusub); ibestUML = ioffs(isub); RAbestUML = RAoffs(rasub);
databestUML = dataout(wsub,nusub,isub,rasub);
UMLbestUML = UMLout(wsub,nusub,isub,rasub);
minoffbestUML = minoffout(wsub,nusub,isub,rasub);
mintbestUML = mintout(wsub,nusub,isub,rasub);
%Print Outputs
fprintf('Max UML Result: %.0e line-km\nUML: %.2f Orbits\nMin. Offset: %.2f m\nTime of Min.: %.0f sec\nwoff: %.3e rad.\nnuoff: %.3e rad.\nioff: %.3e rad.\nRAoff: %.3e rad.\n\n',...
    databestUML,UMLbestUML,minoffbestUML,mintbestdata,wbestUML,nubestUML,ibestUML,RAbestUML)