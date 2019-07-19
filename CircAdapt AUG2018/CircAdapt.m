function CircAdapt
% function CircAdapt
% Core of the CircAdapt model
% Simulates a series of beats based on de contents of structure P
% Theo Arts, Maastricht University, April 26, 2014

global P
P.General.tStart= P.t(end); % time reference for array counter
P.General.tEnd  = P.t(end)+P.General.DtSimulation;
LinkMatrix; % matrices defining connections between elements

% Initialization P.SVar (= state variables) to last state
P2SVar;
P.SVar= P.SVar(end,:);
nSVar = size(P.SVar,2); % start condition
P.t   = P.t(end);
Dt    = P.General.Dt;
% clearing In-Out record, used for Fast adaptation
P.General.In     = [];
P.General.Out    = []; % reset of storage input/output per beat
P.General.FacpControl= 1; % Initialization pCurrent/pTarget
% if FacpControl>1 => loss of blood volume in ArtVen's

Indexation; % Indexation to define connections between elements

ntMax        = ceil((P.General.DtSimulation+P.General.tCycle)/Dt);
SVar         = zeros(ntMax,nSVar); % preallocation SVar
iSVar        = 1; % time counter of SVar
SVar(iSVar,:)= P.SVar; %initial condition
%==== simulation of successive cardiac cycles
while ( P.General.tEnd - P.t(end) > Dt ); 

    % setting time points for a cycle with integer number of sample points
    % by adjustment of cycle length
    nt        = ceil(P.General.tCycle/Dt); % number of time points
    tCycle    = nt*Dt; % set integer number of time steps per cycle
    P.General.tCycle= tCycle; % Adjustment of tCycle to integer time steps
    TimingInit; % Setting depolarization sequence in P.Net.Depolarization
    
    % Allocation additional storage space for delayed tube signals
    itMax  = round((P.t(end)-P.General.tStart+tCycle)/Dt)+1;
    ntExtra= ceil(itMax-size(P.Tube.pP,1)); % Extra needed time points
    % for the upcoming beat    
    Mat    = zeros(ntExtra,P.Tube.n);
    P.Tube.pP = [P.Tube.pP ; Mat]; % proximal source pressure
    P.Tube.pD = [P.Tube.pD ; Mat]; % distal source pressure
    P.Tube.q  = [P.Tube.q  ; Mat]; % mean tube flow
    % existing solution
    disp(['t= ',num2str(P.t(end)),';  Time to go= ',...
        num2str(P.General.tEnd-P.t(end))]); pause(0.01);
    % Differential equation
    SVarAppend    = OdeCA('SVarDot',Dt,tCycle,P.SVar(end,:));
    % Solver ODE especially designed for CircAdapt
    nSVarAppend   = size(SVarAppend,1);
    RgSVar        = iSVar-1+(1:nSVarAppend);
    SVar(RgSVar,:)= SVarAppend; % copy SVarAppend to SVar
    iSVar         = iSVar-1+nSVarAppend; %Last row of SVar
    % SVar(end,:) removed because of overlap with next beat
    P.SVar        = SVarAppend;
    CircDisplay(0); % display and time course of 1-beat hemodynamics
    % Calculation of changes due to adaptation
    Adapt; % Execute Adapt with P.Adapt.FunctionName
    % Time courses in Par belong to parameter setting before Adapt-action!
end

P.SVar= SVar(1:iSVar,:); % State variables of all beats to show trends.
% Be careful, because AdaptXX has changed the P-parameters during the
% simulation so that they do not belong to the Par.SVar.
% Only for Adapt0P (= no adaptation), complete SVar
% is compatible with parameter settings.
end

