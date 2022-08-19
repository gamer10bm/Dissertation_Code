% JULIAN DATE FROM UT DATING SYSTEM
% By: Brian Kaplinger, Ph.D.
% Changed: January 29, 2015
function [jd, j0, ut] = julian_date(y, m, d, hr, min, s)
  % Input Variables:
  %   y: year   (1901 <= y <= 2099)
  %   m: month  (1 <= m <= 12)
  %   d: day    (1 <= d <= 31)
  %   hr: hour  (0 <= hr < 24)
  %   min: minutes (0 <= min < 60)
  %   s: seconds  (0 <= s < 60)
  % Output Variables:
  %   jd:  Julian Date
  %   j0:  Julian Date at 0 hr UT for current date
  %   ut:  Universal Time value in hrs ( < 24.0)
  j0 = 367*y-floor(7*(y+floor((m+9)/12))/4)+floor(275*m/9)+d + 1721013.5;
  ut = hr + min/60 + s/3600;
  jd = j0 + ut/24;
end