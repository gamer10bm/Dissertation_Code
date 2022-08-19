function runjobs(job,waittic,timebetweenupdates)
%Run batch job
    submit(job);
    tic
    lasttoc = toc;
    %Monitor job
    while ~strcmp(job.State,'finished')
        %Hold for a few seconds
        pause(waittic)
        if toc-lasttoc < timebetweenupdates
            lasttoc = toc;
            %Send an update
            numrun = 0;
            numpend = 0;
            numfin = 0;
            for tid = 1:length(job.Tasks)
                %Check for state
                switch job.Tasks(tid).State
                    case 'running'
                        numrun = numrun+1;
                    case 'finished'
                        numfin = numfin+1;
                    case 'pending'
                        numpend = numpend+1;
                end
            end
            
                    
            %Send message
            fprintf('Monte Carlo Job Running:\n\tPending/Finished/Running: %.0f/%.0f/%.0f\n',...
                numpend,numfin,numrun);
        end
    end
end