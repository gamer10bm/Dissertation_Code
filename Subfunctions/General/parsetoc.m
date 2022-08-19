function [elapstr,timeelap, elapunitstr] = parsetoc(t_secs)
%Get time elapsed
timeelap = t_secs;
if timeelap>=60
    %Do hour check
    if timeelap/60>=60
        %Hour
        timeelap = timeelap/60/60;
        elapunitstr = 'hrs';
    else
        %Minute
        timeelap = timeelap/60;
        elapunitstr = 'min.';
    end
else
    %Second
    elapunitstr = 'sec';
end
%Generate time elapsed string
elapstr = sprintf('Time Elapsed: %.2f %s',timeelap,elapunitstr);