function epoch = initial_epoch(yr, month, day, hr, mins, sec)
%  This function stores initial epoch data and initializes site data based
%  on the epoch values given in date and UT.
%  Outputs include julian date, UT in hrs, and sidereal time for Greenwich, 
%  as well as the inputs packaged in case of later use.
%  version A.0; saved B. Kaplinger, 1/4/2020

epoch.year = yr; epoch.month = month; epoch.day = day; epoch.hr = hr;
epoch.min = mins; epoch.sec = sec;

% initial epoch julian date
[jd, j0, ut] = julian_date(yr, month, day, hr, mins, sec);
epoch.JD = jd; epoch.J0 = j0; epoch.UT = ut;

tt0 = (jd-2451545)/36525;
thg0 = (100.4606184+36000.77004*tt0+0.000387933*tt0^2 - ...
    2.583e-8*tt0^3)*pi/180;
thg0 = mod(thg0,2*pi);
epoch.THG0 = thg0;

end