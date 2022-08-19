function [XRK,YRK,ErrorX]=RungeKutta(dervfunc,Ystart,StepSizeorXVec,StartX,EndX,Xexpect)
%RUNGEKUTTA Runs the Runge-Kutta numerical method for scalars and vectors
%   Make sure that the derivative function accepts inputs (t,x) in that
%   order
%   Xexpect is used to estimate the error in the system

%Setup the method call
RKMeth = @(f,h,xn,yn,v1,v2)RKStepVec(f,h,xn,yn);

if nargin<6 %Check if an error needs to be estimated
    Xexpect = 0;
end

if nargin < 4 %Check if a vector is used for X
    %Assumed to be using a vector for X
    Xvec = StepSizeorXVec;
    %Call the method for vectors
    [XRK,YRK]=OneStepMethodXVec(RKMeth,dervfunc,Ystart,Xvec);
    ErrorX = 0;
else %Run the normal method
    [XRK,YRK,ErrorX]=OneStepMethodVec(RKMeth,dervfunc,Ystart,StepSizeorXVec,StartX,EndX,Xexpect);
end
if length(YRK) == 1
    XRK = XRK{1};
    YRK = YRK{1};
end

end