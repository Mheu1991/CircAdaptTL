function TubeDelays
% function StoreDelayedSignals
% Stores sampled pressures for calculation of delayed
% pressures. Needed for Tube wave propagation. Data are collected during
% execution of Ode2 procedure.
% Theo Arts, Maastricht University, Oct 13, 2012

% Imposed changes in the calculation of vicous friction of blood flow, now based on an
% approximation that is in agreement with Womersley theory (considered gold
% standard), Maarten Heusinkveld, August 2018

global P

if P.Tube.n>0 && size(P.t,1)==1; % store signals to be delayed
    Dt=P.General.Dt;
    it=round((P.t-P.General.tStart)/Dt) + 1; % array of time indices
    [nt,nTube]= size(P.Tube.pL);
    DiL       = P.Tube.DiL; % row of delay time intervals L-wave
    DiR       = P.Tube.DiR; % row of delay time intervals R-wave
    itM3      = Mod1(it+(-1:1)',nt);% cyclic time index it+[-1,0,+1]
    itM       = itM3(2,:);          % cyclic time index it

    cL = P.Tube.cL ; % wave velocity to left (backward)
    cR = P.Tube.cR ; % wave velocity to right (forward)
    Len= P.Tube.Len;
    AuxL= Dt*DiL.*cL./Len;
    AuxR= Dt*DiR.*cR./Len;
    DiL = DiL + ( 1 - AuxL ); %delay/Dt left wave
    DiR = DiR + ( 1 - AuxR ); %delay/Dt right wave
    DiL= max(1,min(nt-2,DiL));% safety to extreme delay
    DiR= max(1,min(nt-2,DiR));% safety to extreme delay
    P.Tube.DiL = DiL; % store new delay
    P.Tube.DiR = DiR; % store new delay
    iL= floor(DiL); rL=DiL-iL;
    iR= floor(DiR); rR=DiR-iR;
    FacL= exp(-Dt*DiL.*cL.*P.Tube.Att); % wave attenuation [Att] = 1/m , [cL] = m/s, [Dt] = s , [DiL]= -
    FacR= exp(-Dt*DiR.*cR.*P.Tube.Att); % wave attenuation     

    % Initiation of pressures wave, relative to external tube pressure
    % {pL,pR}= after delay: source pressures {Prox,Dist}
    pTrans= P.Tube.pTrans;
    pM    = P.Tube.p;
    pExt  = pM - pTrans;
    qNP   = P.Tube.qProx; % proximal inflow
    qND   = P.Tube.qDist; % distal outflow
    ZL    = P.Tube.ZL;
    ZR    = P.Tube.ZR;
    pSP   = P.Tube.pProx;
    pSD   = P.Tube.pDist;
    qM    = P.Tube.q(it,:);
    
    ZSum=ZL+ZR;
    pL=-pExt + pSD - qND.*ZSum;
    pR=-pExt + pSP + qNP.*ZSum;
    
    % Flow
    Rgt=it+(0:1);
    dq= 4*((-pL+pR)./ZSum - qM)./(DiL+DiR);
    P.Tube.q(Rgt,:)= [1;1]*(qM+dq);
    % q=flow, ~ synchronous with qProx+qDist
    
    P.Tube.pL(itM3,:)= bsxfun(@times,P.Tube.pL(itM3,:),[0.67;0.5;0])...
        + [0.33;0.5;1]*pL; % itM+1 needed for
    P.Tube.pR(itM3,:)= bsxfun(@times,P.Tube.pR(itM3,:),[0.67;0.5;0])...
        + [0.33;0.5;1]*pR; % itM+1 needed for
    
    pP0= pTrans - qM.*P.Tube.ZR; % DC-correction with attenuation
    pD0= pTrans + qM.*P.Tube.ZL; % DC-correction with attenuation
    % Find indices in cyclic pL and pR memory around interpolation points
    RgL1=Mod1(bsxfun(@plus,(-1:1)',itM-iL),nt);
    RgR1=Mod1(bsxfun(@plus,(-1:1)',itM-iR),nt);
    % Conversion subscript->index
    RgL=bsxfun(@plus,RgL1,(0:nTube-1)*nt);
    RgR=bsxfun(@plus,RgR1,(0:nTube-1)*nt);
    
    % fractional interpolation source pressure
    pLRg=P.Tube.pL(RgL);
    pRRg=P.Tube.pR(RgR);
    pP=pLRg(2:3,:)+bsxfun(@times,conv2(pLRg,[-1;+1],'valid'),rL);
    pD=pRRg(2:3,:)+bsxfun(@times,conv2(pRRg,[-1;+1],'valid'),rR);  
    
    % storage in pP and pD source pressure
    pP02=[pP0;pP0];
    pD02=[pD0;pD0];
    P.Tube.pP(Rgt,:)= [FacL;FacL].*(pP-pP02)+pP02;% correct
    P.Tube.pD(Rgt,:)= [FacR;FacR].*(pD-pD02)+pD02;% approximation    
end

end

function j=Mod1(i,nt)
j=mod(i-1,nt)+1;
end

