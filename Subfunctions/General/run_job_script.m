% script run_job_script
%Send message
elapstr = parsetoc(toc);
fprintf('Submitting Job %.0f of %.0f\n\t%s\n',job_n,totaljobs,elapstr)
%Run batch job
errorsvec = [];
submit(job);
lasttoc = toc;
%Monitor job
printbool = true;
while ~strcmp(job.State,'finished')
    %Hold for a few seconds
    pause(waittic)
    if toc-lasttoc > timebetweenupdates || printbool
        printbool = false;
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
            %Check for errors
            if job.Tasks(tid).hasError && ~any(errorsvec==tid)
                errorsvec(end+1) = tid;
            end
        end
        elapstr = parsetoc(toc);

        %Send message
        fprintf('Job %.0f of %.0f Running:\n\t%s\n\tPending/Finished/Running: %.0f/%.0f/%.0f\n',...
            job_n,totaljobs,elapstr,numpend,numfin,numrun);
        if ~isempty(errorsvec)
            fprintf('\t\t%.0f Errors So Far\n',length(errorsvec))
        end
    end
end
%Ensure everything is done before moving on
wait(job)
%Print success message
elapstr = parsetoc(toc);
fprintf('%s\nJob %.0f of %.0f Completed:\n\t%s\n\t%.0f of %.0f runs completed\n%s\n\n',...
    cmdline,job_n,totaljobs,elapstr,total_n,totalruns,cmdline);
%Check for errors
jobdelvec = [];
for tid = 1:length(job.Tasks)
    if ~job.Tasks(tid).hasError
        %Delete task
        jobdelvec(end+1) = tid;
    end
end
%Delete the tasks that didn't error out
job.Tasks(jobdelvec).delete;
%Send to new job
joberrors = job.recreate;
if ~isempty(job.Tasks)
    %Recreate and resubmit job in case of errors
    fprintf('Resubmitting %.0f task(s)\n',length(job.Tasks))
    submit(joberrors)
    wait(joberrors)
    %Check for errors
    jobdelvec = [];
    for tid = 1:length(job.Tasks)
        if ~joberrors.Tasks(tid).hasError
            %Delete task
            jobdelvec(end+1) = tid;
        end
    end
    %Delete the tasks that didn't error out
    joberrors.Tasks(jobdelvec).delete;
    if ~isempty(joberrors.Tasks)
        error('Check the last job and determine errors before moving forward')
    end
end
%Delete job
delete(job); delete(joberrors);