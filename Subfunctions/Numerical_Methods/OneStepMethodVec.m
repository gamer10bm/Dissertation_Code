function [X,Y,Error,MethodCost] = OneStepMethodVec(MethodStepHandle,FunctionHandle,...
    InitialY,Stepsize,InitialX,GoalX,ExpectedValue)
%This finds y_n+1 using the forward Euler method
    %The function handle should have at least one argument being either x
    %or y
        %If both, use form @(x,y)

dy = FunctionHandle;
%CurrentY is the Y_n at this point or the initial y
y0 = InitialY;
x0 = InitialX;
for hn = 1:length(Stepsize)
    h = Stepsize(hn);
    %Initialize loading vector
    n = 1;
    %Load the first values of X and Y
    Xloop(n) = x0;
    Yloop(:,n) = y0;
    %Start with the n values of x and y
    xn = x0;
    yn = y0;
    %Initialize Cost
    Cost = 0;
    while xn<=GoalX-h+1e-15
        [y_n1,Cloc] = MethodStepHandle(dy,h,xn,yn);
        Cost = Cost+Cloc;
        %Increase loop by 1
        n = n+1;
        yn = y_n1;
        xn = xn+h;
        Yloop(:,n) = yn;
        Xloop(n) = xn;
    end
    for e = 1:numel(ExpectedValue)
        if isnan(ExpectedValue(e))
            Error(e,hn) = nan;
        else
            Error(e,hn) = abs(Yloop(e,end)-ExpectedValue(e));
        end
    end 
    MethodCost(hn) = Cost;
    X{hn} = Xloop;
    Y{hn} = Yloop;
    clearvars Xloop Yloop
end

