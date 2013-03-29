%% Get restoring force for the element
%  returns (6x1) restoring force vector of element in global coordinates
function [Element,RestoringForce] = GetRestoringForce_Type3(Element, Structure, Integrator)

%% Get the displacement Ue (6x1) vector of the element in global coordinates
Ue = zeros(6,1);
dof=zeros(6,1);
dof(1,1)=Element.Nodes(1).UX;
dof(2,1)=Element.Nodes(1).UY;
dof(3,1)=Element.Nodes(1).THETA;
dof(4,1)=Element.Nodes(2).UX;
dof(5,1)=Element.Nodes(2).UY;
dof(6,1)=Element.Nodes(2).THETA;
for i=1:6
    if dof(i,1) ~= -1
        Ue(i,1)=Integrator.Displacement(dof(i,1),1);
    end
end

if( Structure.RigidLinkNo )
    TR = eye(6,6); % Transformation matrix for rigid link
    for jjj=1:2
        idx=find(Structure.RigidLinkNodeID(:,2) == Element.Nodes(jjj).ID); % find if a slave node exists in the element
        if ( idx )
            TL=CalculateRigidLinkTransformation(Structure.RigidLinkMaster{idx},Structure.RigidLinkSlave{idx} );
            TR(3*(jjj-1)+1:3*jjj,3*(jjj-1)+1:3*jjj) = TL;
        end
    end     
    Ue = TR * Ue;          
end

% Input of Drain State Determination Type 2 Element
ndof=6;  % The element has 6 dofs
ddise=Ue-Element.Uprev;  % The incremental displacement vector ddise of the element in global coordinates 
Element.Uprev=Ue; %Update Uprev of element for use in next state determination call

% Calculate element deformation increments
[dvax,dvr] = ELDFO2(ddise,Element.fl,Element.cosa,Element.sina);         

% UPDATE AXIAL FORCES AND AXIAL ENERGY
dfax=Element.eal*dvax;      
Element.ftot(1,1)=Element.ftot(1,1)-dfax;
Element.ftot(2,1)=Element.ftot(2,1)+dfax;
Element.pkg=(Element.ftot(2,1)-Element.ftot(1,1))*0.5;

% MOMENTS IN ELASTIC COMPONENT AT BEGINING OF STEP
bmel(1,1)=Element.bmtot(1,1)-Element.bmep(1,1);
bmel(2,1)=Element.bmtot(2,1)-Element.bmep(2,1);

% LINEAR MOMENT INCREMENTS FOR ELASTO-PLASTIC COMPONENT - SET IN DBM
[dbm] = BMCAO2(dvr,Element); 

% TRACE NONLINEAR PATH FOR ELASTO-PLASTIC COMPONENT---state determination
facac=0.0;
while (facac<0.999999)
    factor=1.0-facac;
    % PLASTIC HINGE ROTATIONS
    kod=Element.kody(1,1)+Element.kody(2,1)-1.0;
    if(kod>0)
        % hinges at both ends
        dpr(1,1)=dvr(1,1);
        dpr(2,1)=dvr(2,1);
    elseif(kod==0)
        % hinge at one end only
        dpr(1,1)=dvr(1,1)+Element.pr21*dvr(2,1);
        dpr(2,1)=dvr(2,1)+Element.pr12*dvr(1,1);
    else
        %never executed
    end
         
    % CALCULATE WHICH END HAS NEXT EVENT IF YIELDS, SET KFAC=IEND AND BBMY=YIELD VALUE
    % SET EVENT FACTOR IN FACTOR
    kfac=3;
    for iend=1:2
        if (Element.kody(iend,1)==0.0)
            dbmm=dbm(iend,1);
            % elastic, get factor for status change(no hinge at end)
            if (dbmm<0.0)
                dbmo=Element.bmy(iend,2)-Element.bmep(iend,1);
                if(dbmo>=dbmm) 
                    fac=dbmo/dbmm;
                else
                    fac=2.0;
                end 
                if (fac<=factor)
                    factor=fac;
                    bbmy=Element.bmy(iend,2);
                    kfac=iend;
                end
            elseif(dbmm>0.0)
                dbmo=Element.bmy(iend,1)-Element.bmep(iend,1);
                if(dbmo<=dbmm) 
                    fac=dbmo/dbmm;
                else
                    fac=2.0;
                end 
                if (fac<=factor)
                    factor=fac;
                    bbmy=Element.bmy(iend,1);
                    kfac=iend;
                end
            else
            end
        else
            % hinge at this end
            if(Element.bmep(iend,1)*dpr(iend,1)<0.0)
                % unloading (plastic work negative)
                factor=0.0;
                kfac=0.0;
                Element.kody(iend,1)=0.0;
            else
            end
        end
    end
    
    % UPDATE BENDING STATE AND CALCULATE BENDING ENERGY
    % dbm= elastic increment - use if elastic bbmy= yield value - use if yields
    if(kfac~=0)
        for iend=1:2
            if(iend==kfac)
                % new hinge at this end
                Element.bmep(iend,1)=bbmy;
                Element.kody(iend,1)=1;
            elseif(Element.kody(iend,1)==0)
                % this end remains elastic
                dmep=factor*dbm(iend,1);
                Element.bmep(iend,1)=Element.bmep(iend,1)+dmep;
            else
                % this end remains plastic
                dppr=factor*dpr(iend,1);
                Element.prtot(iend,1)=Element.prtot(iend,1)+dppr;
                if (dppr<0.0)
                    Element.pracn(iend,1)=Element.pracn(iend,1)+dppr;
                else
                    Element.pracp(iend,1)=Element.pracp(iend,1)+dppr;
                end
            end
        end
    else
    end
        
    % CHECK COMPLETION OF CYCLE
    facac=facac+factor;
    if (facac<0.99999)
        % cycle not completed, calculate new moment increments
        [dbm] = BMCAO2(dvr,Element); 
    end
end

% YIELD MOMENTS FOR NEXT STEP
if(Element.ksurf==1)
    Element = YMOM02(Element);
    % CHECK YIELD SURFACE OVERSHOOT
    for iend=1:2
        if (Element.kody(iend,1)~=0.0)
            % inelastic on old surface
            if (Element.bmep(iend,1)>Element.bmy(iend,1)) 
                % inelastic positive on new surface also
                Element.bmep(iend,1)=Element.bmy(iend,1);
            elseif (Element.bmep(iend,1)<Element.bmy(iend,2))
                % inelastic negative on new surface also
                Element.bmep(iend,1)=Element.bmy(iend,2);
            elseif ((Element.bmep(iend,1)<(Element.bmy(iend,1)*0.98)) && (Element.bmep(iend,1)>(Element.bmy(iend,2)*0.98)))
                % moves inside new surface by > 2% - reset yield code
            Element.kody(iend,1)=0;
            else
            end
        else
            % elastic on old surface
            if (Element.bmep(iend,1)>Element.bmy(iend,1))
                % inelastic positive on new surface
                Element.bmep(iend,1)=Element.bmy(iend,1);
                Element.kody(iend,1)=1;
            elseif (Element.bmep(iend,1)<=Element.bmy(iend,2))
                % inelastic negative on new surface
                Element.bmep(iend,1)=Element.bmy(iend,2);
                Element.kody(iend,1)=1;
            else
            end
        end
    end
end

% UPDATE TOTAL FORCES
dsf=(Element.bmtot(1,1)+ Element.bmtot(2,1))/Element.fl;
ppsh=Element.psh/(1.0-Element.psh);
dmeli=(Element.ek11*dvr(1,1)+Element.ek12*dvr(2,1))*ppsh+bmel(1,1);
dmelj=(Element.ek12*dvr(1,1)+ Element.ek22*dvr(2,1))*ppsh+bmel(2,1);
Element.bmtot(1,1)=Element.bmep(1,1)+dmeli;
Element.bmtot(2,1)=Element.bmep(2,1)+dmelj;
dsf=dsf-(Element.bmtot(1,1)+Element.bmtot(2,1))/Element.fl;
Element.sftot(1,1)=Element.sftot(1,1)-dsf;
Element.sftot(2,1)=Element.sftot(2,1)+dsf;
RestoringForce = ELASO2(ndof,Element);

if( Structure.RigidLinkNo )
    RestoringForce = TR' * RestoringForce;  
end


function [dbm] = BMCAO2(dvr,vars)
dbm=zeros(2,1);
% SET KYY TO INDICATE TYPE OF ELEMENT STIFFNESS TO USE
kyy=vars.kody(1,1)+2*vars.kody(2,1)+1;
if(kyy==1)
    % No hinge at either end
    dbm(1,1)=vars.ek11*dvr(1,1)+vars.ek12*dvr(2,1);
    dbm(2,1)=vars.ek12*dvr(1,1)+vars.ek22*dvr(2,1);
elseif(kyy==2)
    % Hinge at end 1
    dbm(1,1)=0.0;
    dbm(2,1)=vars.ek22h*dvr(2,1);
elseif(kyy==3)
    % Hinge at end 2
    dbm(1,1)=vars.ek11h*dvr(1,1);
    dbm(2,1)=0.0;
elseif(kyy==4)
    % Hinges at both end
    dbm(1,1)=0.0;
    dbm(2,1)=0.0;
else
end
	  	 
function [relas] = ELASO2(ndof,vars)
relas=zeros(ndof,1);
% BASIC RESISTING FORCES
relas(1,1)=(vars.ftot(1,1)-vars.aint(1,1))*vars.cosa-(vars.sftot(1,1)-vars.sfint(1,1))*vars.sina;
relas(2,1)=(vars.ftot(1,1)-vars.aint(1,1))*vars.sina+(vars.sftot(1,1)-vars.sfint(1,1))*vars.cosa;
relas(3,1)=vars.bmtot(1,1)-vars.bmint(1,1)-vars.beint(1,1);
relas(4,1)=(vars.ftot(2,1)-vars.aint(2,1))*vars.cosa-(vars.sftot(2,1)-vars.sfint(2,1))*vars.sina;
relas(5,1)=(vars.ftot(2,1)-vars.aint(2,1))*vars.sina+(vars.sftot(2,1)-vars.sfint(2,1))*vars.cosa;
relas(6,1)=vars.bmtot(2,1)-vars.bmint(2,1)-vars.beint(2,1);

function [dvax,dvr] = ELDFO2(ddise,f1,cosa, sina) 
dvr=zeros(2,1);
% DEFORMATION INCREMENTS OR DEFORMATION RATES
dvax=cosa*(ddise(4,1)-ddise(1,1))+sina*(ddise(5,1)-ddise(2,1));
rot=(sina*(ddise(4,1)-ddise(1,1))+cosa*(ddise(2,1)-ddise(5,1)))/f1;
dvr(1,1)=ddise(3,1)+rot;
dvr(2,1)=ddise(6,1)+rot;

function vars = YMOM02(vars)
facc=-1.0;
for iend = 1 : 2
    % Change sign of axial force and end 1
    fft = vars.ftot(iend,1) * facc;
    facc = 1;
    fac = 1;
    for j = 1 : 2
        fac = -fac;
        yieldcalc = 0;  % Flag if the yield moment is calculated inside 
                        % this next for loop
        for i = 1 : 3
            if (fft < vars.pmx(i,j,iend))
                yieldcalc = 1;  % Set flag since the condition was met
                % Calculate yield moment (M = a2 + a1P)
                bbmy = vars.a2(i,j,iend) - vars.a1(i,j,iend) * fft;
                if (fac * bbmy > 0) 
                    % Axial force is beyond Pc and Pt
                    bbmy = 0;
                end
                vars.bmy(iend,j) = bbmy;
            end
        end
        % Only calculate yield moment if it wasn't calculated in the above
        % for loop
        if (yieldcalc == 0)
            i = 4;
            % Calculate yield moment (M = a2 + a1P)
            bbmy = vars.a2(i,j,iend) - vars.a1(i,j,iend) * fft;
            if (fac * bbmy > 0) 
                % Axial force is beyond Pc and Pt
                bbmy = 0;
            end
            vars.bmy(iend,j) = bbmy;
        end
    end
end
fac=-vars.bmy(1,1);
vars.bmy(1,1)=-vars.bmy(1,2);
vars.bmy(1,2)=fac;