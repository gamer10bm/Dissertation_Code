%% Settings
%Radar Settings 
rad_freq = 147e6; %Hz
clight = physconst('Lightspeed');

%Valid Array Spacing
avg_array_minwl = 0.1; %wavelengths
avg_array_maxwl = 2.5; %wavelengths
avg_array_rangewl = [avg_array_minwl avg_array_maxwl];

%% Automated Calculations
wavelength = clight/rad_freq; %m