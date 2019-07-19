function MemAlloc
% function MemAlloc
% Allocates memory for matrices with column length= number of time points
% These allocations are necessary because elements in matrices get a value
% by gradual substitutions, using indices for locating them in the matrix.
% Theo Arts, Maastricht University, Oct 13, 2012

global P;

nt=size(P.t,1); % number of time points

% test on change of column length, only then MemAlloc is active
ExecuteMemAlloc=1;
if isfield(P,'SVarDot')
    ExecuteMemAlloc = (size(P.SVar,1)~=size(P.SVarDot,1));
end
    
if ExecuteMemAlloc
    
    nWall  =P.Wall.n; % number of walls
    nCavity=P.Cavity.n; % number of cavities
    nPatch =P.Patch.n; % number of patches
    nNode  =P.Node.n; % number of nodes
    nTriSeg=P.TriSeg.n; % number of TriSeg objects
    nBag   =P.Bag.n; % number of TriSeg objects
    
    % Matrix constructions
    MatWall  = zeros(nt,nWall);
    MatCavity= zeros(nt,nCavity);
    MatPatch = zeros(nt,nPatch);
    MatNode  = zeros(nt,nNode);
    MatTriSeg= zeros(nt,nTriSeg);
    MatBag   = zeros(nt,nBag);
    
    % Allocating matrices for the walls
    P.Wall.Am0   =MatWall; %unstressed mid-wall area
    P.Wall.DADT  =MatWall; %area compliance
    P.Wall.T     =MatWall; %tension
    P.Wall.Cm    =MatWall; %curvature
    P.Wall.Am    =MatWall; %actual mid-wall area
    P.Wall.pTrans=MatWall; %transmural pressure

    % Allocating matrices for the cavities
    P.Cavity.A     =MatCavity; %cross-sectional area
    P.Cavity.Z     =MatCavity; %resistive impedance
    P.Cavity.pTrans=MatCavity; %resistive impedance

    % Allocating matrices for the patches
    P.Patch.T=MatPatch; %Tension

    % Allocating matrices for the nodes
    P.Node.q=MatNode; %node flow source
    P.Node.Y=MatNode; %node internal flow conductance

    % Allocating matrices for the TriSeg objects
    P.TriSeg.VS  =MatTriSeg; % volume displacement septal wall
    P.TriSeg.YS  =MatTriSeg; % radius junction ring
    P.TriSeg.VDot=MatTriSeg; % SVar-Dot
    P.TriSeg.YDot=MatTriSeg; % SVar-Dot
    
    P.Bag.V = MatBag;
end
end

