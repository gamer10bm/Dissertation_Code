%% Initialize STK
try
    %Grab open Instance
    app = actxGetRunningServer('STK12.application');
catch
    %Grab STK Instance
    app = actxserver('STK12.application');
end
%Change control and visibility
usercontrol = true;
app.UserControl = usercontrol;
app.Visible = usercontrol;
%Grab application root
root = app.Personality2;

%Check if scenario proper scenario is open
scencheck = root.Children.Contains('eScenario',scen_name);
if scencheck
    %Get current Scenario
    scenario = root.Children.Item(scen_name);
    %Unload satellite objects
    
else
    %Close any open scenario
    root.CloseScenario
    %Check if scenario already exists
    existscen = exist(savepath_stk,'file');
    if 0 || ~existscen
        %Create Scenario
        scenario = root.Children.New('eScenario',scen_name);
    else
        %Open Scenario (debug)
        root.LoadScenario(savepath)
        %Grab the current scenario
        scenario = root.CurrentScenario;
    end
end
%Set the analytical time period
scenario.SetTimePeriod(timestartstr,sprintf('+%0.0fhr',daystoend*24));
%Reset the animation time
root.ExecuteCommand('Animate * Reset');