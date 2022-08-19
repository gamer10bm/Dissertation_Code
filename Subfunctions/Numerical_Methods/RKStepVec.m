function [y_n1,Cost] = RKStepVec(FunctionHandle,Stepsize,CurrentX,CurrentY)
%This finds y_n+1 using the forward Euler method
    %The function handle should have at least one argument being either x
    %or y
        %If both, use form @(x,y)
dy = FunctionHandle;

%Be careful about this one
h = Stepsize;
%CurrentY is the Y_n at this point or the initial y
yn = CurrentY;
xn = CurrentX;
xn_1 = xn+h;
xn_1_2 = xn+.5*h;

%Perform RK calculations
k1 = h.*dy(xn,yn);
k2 = h.*dy(xn_1_2,yn+.5*k1);
k3 = h.*dy(xn_1_2,yn+.5*k2);
k4 = h.*dy(xn_1,yn+k3);
    
y_n1 = yn+(1/6).*(k1+2.*k2+2.*k3+k4);
Cost = 4;
end
