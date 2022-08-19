function [JD_UTI] = UTC2JD(yr,mnth,day,hr,min,sec)
%
% Function based on algorithm 14 of Vallado
%
% Created by: Bailey Miller


JDyr = 367*yr;
JDmonth = -floor((7*(yr+floor((mnth+9)/12)))/4);
JDmonth2 = floor(275*mnth/9);
JDhrminsec = ((sec/60+min)/60+hr)/24;

JD_UTI = JDyr+JDmonth+JDmonth2+day+1721013.5+JDhrminsec;
end