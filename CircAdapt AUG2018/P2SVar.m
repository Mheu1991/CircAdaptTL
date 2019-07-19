function P2SVar
% function P2SVar
% Structure P -> scaling -> State variables stored in P.SVar
% Time in column direction, SVar in row direction
% Theo Arts, Maastricht University, Jan 5, 2014

global P

ScaleVqY=P.General.ScaleVqY;
ScV=ScaleVqY(1); Scq=ScaleVqY(2);

% Conversion of volumes, flows and distances to conveniently scaled
% variables. To avoid scaling, volumes are converted to its logarithm.
% Variables crossing zero (Flows, TriSegV) are converted to its asinh.
LnVTube  = log(P.Tube.V);
LnVCavity= log(P.Cavity.V);
LnTriY   = log(P.TriSeg.Y);
AsinhqValve= asinh(P.Valve.q/Scq);
AsinhTriV  = asinh(P.TriSeg.V/ScV);
P.SVar=[P.t, LnVCavity, LnVTube, AsinhqValve ,...
    P.Patch.C,P.Patch.Lsi, AsinhTriV, LnTriY];
end

