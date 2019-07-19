function FI=OdeCA(DfDtFu,dt,TMax,fp)
% function [FI]=OdeCA(DfDtFu,dt,TMax,fp)
% Solving system of ordinary differential equation with given steps,
% 2nd order, 1 calculation of derivative per time step
% Especially designed for CircAdapt, variable ending of integration
% Thus, changes in heart rate can be simulated
% 1st state variable fi(1) is time with [dfi/dt(1)=DfDtFu(1)]=1
% If DfDtFu(1)==0, sign to end integration before TMax is reached
% DfDtFu= string pointing to function df/dt
% dt= integration step
% TMax= maximum duration of integration
% OUTPUT: t=column vector of time points
% FI: matrix, rows -> time : columns -> state variables
% If integration ends before TMax, FI is shorter in rows
% Theo Arts, April 26, 2014
dfdt   = str2func(DfDtFu); % setting derivative function
nt     = round(TMax/dt)+1; % max number of time points
np     = numel(fp);        % number of state variables
fi     = zeros(nt,np);     % reserve memory space for output
h      = fp;               % h=output, set to initial condition
fi(1,:)= h ;               % copy to output matrix
t      = fp(1);            % First state variable represents time
df1    = dfdt(t,h',[])';   % 1st derivative calculation
df1(1) = 1;
f      = h + df1*dt;
t      = t+dt;
df     = dfdt(t,f',[])';   % 1st derivative calculation
df(1)  = 1;
g      = 0.5*(df+df1);
h      = h + g*dt;         % value of state variable, 2nd order approx
fi(2,:)= h;                % copy to output matrix
it=3;                      % Regular integration starts with 3rd point
while it<=nt               % stops integration
    f = h + g*dt;
    t = t+dt;
    df= dfdt(t,f',[])';    % derivative calculation
    if df(1) == 0          % df(1)==0 is used as sign to stop integration
        nt=it;             % stop integration
        df(1) = 1;
    end
    g = (g+2*df)/3;        % updated derivative
    h = h + g*dt;          % value of state variable, 2nd order approx
    fi(it,:)=h;            % copy to output matrix
    it=it+1;
end
FI=fi(1:nt,:);
end

