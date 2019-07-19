function Saturation
global P;

% Oxygen information
P.General.Oxygen.Hb=15;
P.General.Oxygen.k =1.34;
%Standard O2 content in literature in [mlO2/100ml] blood

% Venous Target concentrations are stored in P.ArtVen.SatVen
nBeat=P.General.nBeatSat; % number of cycles for to approach steady state of saturation

% Calculation of Node saturations
P.Node.V=P.Cavity.V*P.Cavity.Cavity2Node+...
    P.Tube.V*(0.5*(P.Tube.Tube2NodeProx+P.Tube.Tube2NodeDist));
% volumes attributed to nodes

NodeSat=ones(P.Node.n,1); % row vector of initial saturation in nodes
if isfield(P.Node,'Sat') && size(P.Node.Sat,2)== P.Node.n
    NodeSat= P.Node.Sat(end,:)'; % row vector of initial saturation in nodes
end
nt=length(P.t); % number of samples in time
Aux=zeros(size(P.Node.V)); % matrix of time-dependent node variables
P.Node.Sat   = Aux; % saturation-volume product ~ oxygen content
P.Node.SatDot= Aux; % saturation-volume product ~ oxygen content

S=zeros(nBeat*nt,P.Node.n); % used for display of saturation in nodes
for i=1:nBeat
    it= 1:nt;
    [tDummy,Sat]= Ode2('SatDit',...
        it,NodeSat); % solving of Differential Equations
    SatDit(1:nt,Sat',[]); % substitution of solution
    S((i-1)*nt+(1:nt),:)=P.Node.Sat; % filling display matrix beat by beat
    NodeSat=P.Node.Sat(end,:)'; % new starting condition
end
t=P.General.tCycle*(1:nBeat*nt)'/nt; % time vector of samples for display
figure(2); plot(t,S) % plot of saturation as fu of time
% if mean(imag(S(:)))~=0
%     disp(stop)
% end

end

%=============== Odesolver and Derivative  =========================

function [ti,fi]=Ode2(DfDtFu,ti,fp)
% function [ti,fi]=Ode2(DfDtFu,ti,fp)
% Solving system of ordinary differential equation with given steps, 2nd order
% 1 calculation of derivative per time step
% DfDtFu= string pointing to function df/dt
% ti= time points
% fp= row vector of state variables defines initial condition at time ti(1)
nt=numel(ti); % number of time points
np=numel(fp); % number of state variables
fi=zeros(nt,np); % reserve memory space for output
fi(1,:)=fp; % initial condition
dfdt=str2func(DfDtFu); % setting derivative function
f1=fp; % initial condition
df1= dfdt(ti(1),fp',[])'; % 1st derivative calculation
for it=2:nt % step by step integration
    dt=ti(it)-ti(it-1); % time step
    fa=f1+dt*df1; % estimate next point
    df2= dfdt(ti(it),fa',[])'; % derivative calculation
    f2=f1+(dt/2)*(df1+df2); % value of state variable, 2nd order approx
    fi(it,:)=f2; % copy to output matrix
    f1=f2; df1=df2; % copy SVar-value and derivative for next step
end

end

function SatDitT= SatDit(it,NodeSat,flag)
% function SatDitT=SatDit(it,NodeSat,flag)
% Derivative SatDitT of product Sat=Saturation, to index increment
% it      = array index of time points, proportional to time
% NodeSat= column (nodes) of state variables Sat (P.Node.Sat),
% proportional to oxygen content
% flag    = not used
% This function is suited to be solved with Ode2, the solution is
% equidistantly sampled

global P;

%SatVeTarget=P.General.Oxygen.SatVen; % Venous target saturation of blood
Dt   = P.General.Dt; % time increment
nt   = length(it); % number of time points
Col1t= ones(nt,1);

P.Node.Sat(it,:)= NodeSat';

%Valves
qV   = P.Valve.q(it,:);
DSatV= P.Node.Sat(it,P.Valve.iNodeProx) - P.Node.Sat(it,P.Valve.iNodeDist);

%ArtVen elements
qAv        = P.ArtVen.q(it,:);
iAvNodeArt= P.Cavity.iNode(P.ArtVen.iNode);
iAvNodeVen= P.Cavity.iNode(P.ArtVen.iNode+1);
% Oxygen consumption/uptake so that ArtVen outflow concentration == Target
DSatAvFwd= Col1t*P.ArtVen.SatVen - P.Node.Sat(it,iAvNodeVen);
DSatAvBwd= P.Node.Sat(it,iAvNodeArt) - Col1t*P.ArtVen.SatVen;% never happens?

%Tube elements
qT        = 0.5*(P.Tube.qProx(it,:)+P.Tube.qDist(it,:));
iTNodeProx= P.Tube.iNodeProx;
iTNodeDist= P.Tube.iNodeDist;
DSatT     = P.Node.Sat(it,iTNodeProx) - P.Node.Sat(it,iTNodeDist);

%time derivative of oxygen Saturation in nodes
P.Node.SatDot(it,:)= (...
    (max(0,qV).*DSatV)*P.Valve.Valve2NodeDist + ...
    (min(0,qV).*DSatV)*P.Valve.Valve2NodeProx + ...
    (max(0,qAv).*DSatAvFwd)*P.ArtVen.ArtVen2NodeVen + ...
    (min(0,qAv).*DSatAvBwd)*P.ArtVen.ArtVen2NodeArt + ...
    (max(0,qT ).*DSatT )*P.Tube.Tube2NodeDist + ...
    (min(0,qT ).*DSatT )*P.Tube.Tube2NodeProx ...
    )./P.Node.V(it,:);

SatDitT=Dt*P.Node.SatDot(it,:);

end

