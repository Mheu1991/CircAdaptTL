function OutNew=SteadyStateP
% function OutNew=SteadyStateP
% OutNew= new estimate vector of steady state parameters, to be used
% as start condition for next heart beat
% Theo Arts, Maastricht University, Oct 13, 2012

global P

% read input-output record, take logarithm for relative change
In  = log(P.General.In );% values at input
Out = log(P.General.Out);% values after CircAdapt
[nx,np]= size(In); % [number of samples, number of parameters]
yRef= Out(end,:); % last y-value vector is used as reference

nxn = max(1,mod(np,3)); % number of vectors used to predict next solution
nMin= max(2,nx-nxn); % number of samples used for prediction
Rgn =(nMin:nx)'; % indices used for prediction
nxn = length(Rgn);
if nxn>1
    nt = size(Rgn,1); Col1=ones(nt,1);
    X  = bsxfun(@plus,In(Rgn,:) ,-yRef); %input
    Y  = bsxfun(@plus,Out(Rgn,:),-yRef); %output
    YmX= Y-X;
    Ar = double(~cellfun(@isempty,strfind(P.Cavity.Name,'Ar')));
    Ve = double(~cellfun(@isempty,strfind(P.Cavity.Name,'Ve')));
    Pu = double(~cellfun(@isempty,strfind(P.Cavity.Name,'Pu')));
    A  = bsxfun(@times,[Ar;Ve],1-Pu)'; % Systemic arteries and vein
    M1 = [YmX*A,Col1];
    NY0= pinv(M1)*(Y*A);
    dyA= NY0(end,:);
    dy = dyA*pinv(A); % estimated correction
    yNew= yRef+dy; % reasonable estimate
else
    yNew=yRef;
end

OutNew=exp(yNew); % prediction after exponential

end

