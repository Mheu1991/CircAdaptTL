function P2SVarDot
% function P2SVarDot
% Structure P -> scaling -> d/dt State variables stored in P.SVarDot
% Time in column direction, SVar in row direction
% Theo Arts, Maastricht University, Jan 5, 2014

global P

ScaleVqY=P.General.ScaleVqY;
ScV=ScaleVqY(1); Scq=ScaleVqY(2); ScY=ScaleVqY(3); Sc0=1;

% Conversion of volumes, flows and distances to conveniently scaled
% variables. To avoid scaling, volumes are converted to its logarithm.
% Variables crossing zero (Flows, TriSegV) are converted to its asinh.

LnVTubeDot  = P.Tube.VDot   ./ P.Tube.V  ;
LnVCavityDot= P.Cavity.VDot ./ P.Cavity.V;
LnTriYDot   = P.TriSeg.YDot ./ P.TriSeg.Y;
AsinhqValveDot= P.Valve.qDot /Scq ./ sqrt(1+(P.Valve.q/Scq ).^2);
AsinhTriVDot  = P.TriSeg.VDot/ScV ./ sqrt(1+(P.TriSeg.V/ScV).^2);
P.SVarDot= [P.tDot,LnVCavityDot,LnVTubeDot, ...
    AsinhqValveDot, P.Patch.CDot, P.Patch.LsiDot, ...
    AsinhTriVDot, LnTriYDot];

end

