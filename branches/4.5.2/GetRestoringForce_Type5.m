%% Get restoring force for the element
%  returns (6x1) restoring force vector of element in global coordinates
% Written by: Theodore L. Karavasilis 
% returns the resisting force vector and also calculates the current tangent
% stiffness 
function [Element,RestoringForce] = GetRestoringForce_Type5(Element,Structure,Integrator)
% Get the displacement Ue (6x1) vector of the element in global coordinates
Uev = zeros(6,1);
dof=zeros(6,1);
dof(1,1)=Element.Nodes(1).UX;
dof(2,1)=Element.Nodes(1).UY;
dof(3,1)=Element.Nodes(1).THETA;
dof(4,1)=Element.Nodes(2).UX;
dof(5,1)=Element.Nodes(2).UY;
dof(6,1)=Element.Nodes(2).THETA;
for i=1:6
    if dof(i,1) ~= -1
        Uev(i,1)=Integrator.Displacement(dof(i,1),1);
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
    Uev = TR * Uev;         
end

% Total displacement (the element is used as a rotational spring)-if you
% to use it as translational horizontal, then simply
% Ue=Uev(1,1)-Uev(4,1)......e.t.c
Ue = Uev(3,1)-Uev(6,1); 

% Incremental displacement
ddise=Ue-Element.Uprev;  %the incremental displacement of the element 

% Update Uprev
Element.Uprev=Ue; 

% Calculation of deltaD
deltaD=ddise;
d=Ue; 
Element.dP = d-deltaD;	   
if (d>0.0) 
    [temp_1_farzin]=interPoint(0.0,Element.fCapRefPos,Element.capSlope*Element.elstk,0.0,Element.Resfac*Element.fyieldPos,0.0); 
    if (d<temp_1_farzin) 
        Element.iNoFpos = 0;
        Element.LP=0;
    end
else
    [temp_1_farzin]=interPoint(0.0,Element.fCapRefNeg,Element.capSlopeNeg*Element.elstk,0.0,Element.ResfacNeg*Element.fyieldNeg,0.0); 
    if (d>temp_1_farzin) 
        Element.iNoFneg = 0;
        Element.LN=0;
    end
end

% Other variables
Element.flagdeg = 0;
betas = 0.0;
betad = 0.0;
% Yield displacements
Element.dyieldPos = Element.fyieldPos/Element.elstk;
Element.dyieldNeg = Element.fyieldNeg/Element.elstk;

% Initialize parameters in the first cycle
if (Element.kon==0) 
    ekhard = Element.elstk*Element.alpha;
    Element.alphaNeg=Element.alphaN;
    Element.alphaPos=Element.alpha;
    Element.ekhardPos = ekhard;
    Element.ekhardNeg = ekhard;
    Element.dyieldPos = Element.fyieldPos/Element.elstk;
    Element.dyieldNeg = Element.fyieldNeg/Element.elstk;
    Element.Enrgts = Element.fyieldPos*Element.dyieldPos*Element.ecaps;
    Element.Enrgtk = Element.fyieldPos*Element.dyieldPos*Element.ecapk;
    Element.Enrgtd = Element.fyieldPos*Element.dyieldPos*Element.ecapd;
    Element.dmax = Element.dyieldPos;
    Element.dmin = Element.dyieldNeg;     
    Element.ekunload = Element.elstk;
    Element.ekexcurs = Element.elstk;
    Element.Enrgtot = 0.0;
    Element.Enrgc = 0.0;
    Element.fyPos = Element.fyieldPos;
    Element.fyNeg = Element.fyieldNeg;
    Element.dyNeg=Element.dyieldNeg;
    Element.dyPos=Element.dyieldPos;
    Element.resSn = Element.fyieldPos;
    Element.resSp = Element.fyieldNeg;
    Element.cpPos = Element.capDispPos;
    Element.cpNeg = Element.capDispNeg;
    fPeakPos=Element.fyieldPos+ekhard*(Element.capDispPos-Element.dyieldPos);
    fPeakNeg=Element.fyieldNeg+ekhard*(Element.capDispNeg-Element.dyieldNeg);
    Element.flagControlResponse = Element.flagControlResponse;      % This flag Controls the ultimate rotation 
    Element.DPlus = Element.DPlus;                                  % Slab Effect Positive Direction
    Element.DNeg = Element.DNeg;                                    % Slab Effect Negative Direction
    if (Element.cpPos<Element.dyieldPos) 
        fPeakPos =Element.fyieldPos*Element.cpPos/Element.dyieldPos;
    end
    if (Element.cpNeg>Element.dyieldNeg) 
        fPeakNeg =Element.fyieldNeg*Element.cpNeg/Element.dyieldNeg;
    end
    Element.fCapPos = fPeakPos;
    Element.fCapNeg = fPeakNeg;   
    Element.LP = 0;
    Element.LN = 0;
    Element.fLimPos = 0;
    Element.fLimNeg = 0;
    Element.dLimPos = 0;
    Element.dLimNeg = 0;		
    Element.iNoFpos = 0;
    Element.iNoFneg = 0;
    Element.interup=0;
    Element.fCapRefPos=-Element.capSlope*Element.elstk*Element.capDispPos+fPeakPos;
    Element.fCapRefNeg=-Element.capSlopeNeg*Element.elstk*Element.capDispNeg+fPeakNeg;
    Element.capSlopeOrig = Element.capSlope;
    Element.capSlopeOrigNeg = Element.capSlopeNeg;
    Element.flagstopdeg	= 0;
    Element.dCap2Pos = 0.0;
    Element.dCap1Pos = 0.0;
    Element.dCap2Neg = 0.0;
    Element.dCap1Neg = 0.0;
    Element.dyPos = 0.0;
    Element.dyNeg = 0.0;
    Element.RSE = 0.0;    
    Element.ekt = 0.0;   
    if(deltaD>=0.0)
        Element.kon = 1;
    else
        Element.kon = 2;
    end
end

% Variables set to appease conditional statements for big loop
Element.fyNeg = Element.fyNeg;
Element.fyPos = Element.fyPos;
Element.ekhardNeg = Element.ekhardNeg;
Element.ekhardPos = Element.ekhardPos;
Element.dmin = Element.dmin;
Element.dmax = Element.dmax;
Element.resSn = Element.resSn;
Element.resSp = Element.resSp;
Element.ekunload = Element.ekunload;
Element.alphaPos = Element.alphaPos;
Element.alphaNeg = Element.alphaNeg;
Element.Enrgtk = Element.Enrgtk;
Element.ekexcurs = Element.ekexcurs;
Element.cpPos = Element.cpPos;
Element.cpNeg = Element.cpNeg;

%%
% 	******************* S T A R T S   B I G   L O O P  ****************
% 	IF   D E L T A > 0 - - - - - - - - - - - - - - - - - - - - - - - -  
if (deltaD>=0.0) 
    if (Element.iNoFpos==1) 
        [Element.dNewLoadPos]=interPoint(Element.dyNeg,Element.fyNeg,Element.ekhardNeg,0.0,0.0,0.0);
    end
 	%If there is a corner changing delta from negative to positive      
    if (Element.kon==2)       
		Element.kon = 1;
		Element.dlstNeg = Element.dP;
		Element.flstNeg = Element.fP;

        if (Element.dlstNeg<=Element.dmin)
			Element.fLimNeg = Element.flstNeg;
			Element.dLimNeg = Element.dlstNeg;
        end
     
        if (Element.resSn>0)
			Element.RSE = 0.5*Element.fP*Element.fP/Element.ekunload;
		else
			Element.RSE=0.5*(Element.fP+Element.resSn)*(Element.sn-Element.dP)+0.5*Element.resSn*Element.resSn/(Element.elstk*Element.alphaPos);
        end               
        
        if((Element.ecapk~=0.0) && (d<Element.sp) && abs(Element.capSlope)>=1.0e-3 && (abs(Element.capSlopeNeg)>=1.0e-3))
            betak = ((Element.Enrgc-Element.RSE)/(Element.Enrgtk-(Element.Enrgtot-Element.RSE)))^Element.ck;
            if(((Element.Enrgtot-Element.RSE)>=Element.Enrgtk)||(betak>=1.0)) 
                betak = 1.0;
            end
            if( Element.flagstopdeg ~=1 && Element.flagdeg ~=1) 	        
                Element.ekunload = Element.ekexcurs*(1-betak);
            else
                Element.ekunload = Element.ekexcurs;
            end
            if(Element.ekunload<=0.1*Element.elstk) 
                Element.ekunload = 0.1*Element.elstk;
            end
        end

        % Element.sn CALCULATION-------------------------------------------------
        if((Element.dmin<Element.dyieldNeg)||(Element.dmax>Element.dyieldPos)) 
            if((Element.dP<Element.sp)||(((Element.dP>Element.sp)&&(Element.ekP==Element.ekunload))))
                [Element]=snCalc(Element);

                if((abs(Element.dmax-Element.dyieldPos)>=1.0e-10)&&(abs(Element.sn-Element.dyieldPos)<=1.0e-10)) 
                    Element.sn=Element.dyieldPos-1.0e-9;
                end
            end
        end
       
        if (Element.sn > Element.dmax)
            Element.dmax = Element.sn;
        end
        
        if ((Element.iNoFneg==1)&&(Element.dP<=Element.dNewLoadNeg)) 
			Element.sn = Element.dNewLoadNeg - 1.0e-6;
			Element.resSn = 0;
        end
    end
    % 	LOADING ----------------------------------------------------------
    %      Push envelope
    % Compute force and stiffness for this increment
    if ((Element.iNoFneg==1)&&(Element.iNoFpos==1)&&(d<Element.fracDispPos)) 
        f = Element.Resfac*Element.fyieldPos;
		Element.ek = 0.0;		
    elseif((Element.iNoFneg==1)&&(Element.iNoFpos==1)&&(d>=Element.fracDispPos)...
            ||Element.flagControlResponse==1) 
	    f = 1.0e-10;
	    Element.ek = 0.0;	    
	    Element.flagstopdeg = 1;
	elseif ((Element.iNoFpos==1)&&(d>Element.sn)&&(d<Element.fracDispPos)) 
		f = Element.Resfac*Element.fyieldPos;
		Element.ek = 0.0;		
	elseif ((Element.iNoFpos==1)&&(d>Element.sn)&&(d>=Element.fracDispPos))
		f = 1.0e-10;
		Element.ek =0.0;		
	    Element.flagstopdeg	= 1;
    elseif (d>=Element.fracDispPos || Element.flagControlResponse == 1)
		f = 1.0e-10;
		Element.ek =0.0;		
		Element.flagstopdeg = 1;
        Element.flagControlResponse = 1;                 % Switches the response to zero strength
    elseif ((Element.iNoFpos==1)&&(d<Element.sn)) 
		Element.ek = Element.ekunload;
		f  = Element.fP+Element.ek*deltaD;		
	elseif ((Element.iNoFneg==1)&&(d<Element.dNewLoadNeg))
		f = 0;
		Element.ek = 0;		
	elseif (d>Element.dmax)
        [d,f,Element.ek]=envelPosCap2(Element.fyPos,Element.alphaPos,Element.capSlope,Element.cpPos,d,0.0,Element.ek,Element.elstk,Element.fyieldPos,Element.Resfac);
		Element.dmax = d;				                

        % COMPUTE MAXIMUM POSSIBLE DISPLACEMENT        
        [dBoundPos,Element]=boundPos(Element);
        if((d>dBoundPos)||(Element.ek==1.0e-7)) 
            Element.iNoFpos = 1;
        end
	elseif (abs(Element.sn)>1.0e-10) 
        if (Element.LP==0) 
            if (Element.cpPos<=Element.dyPos)
                [xDevPos1,yDevPos1]=interPoint(Element.sn,Element.resSn,Element.ekhardPos,Element.cpPos,Element.fCapPos,Element.capSlope*Element.elstk);
                [xDevPos2,yDevPos2]=interPoint(Element.sn,Element.resSn,Element.ekunload,Element.cpPos,Element.fCapPos,Element.capSlope*Element.elstk);
                if (xDevPos1>xDevPos2) 
                    xDevPos = xDevPos1;
                    yDevPos = yDevPos1;
                else
					xDevPos = xDevPos2;
					yDevPos = yDevPos2;
                end		

                if ((d<=Element.sn)&&(d<=xDevPos)) 
 					Element.ek = Element.ekunload;
					f  = Element.fP+Element.ek*deltaD;					
                elseif ((d>Element.sn)&&(d<=xDevPos)) 
					Element.ek = Element.elstk*Element.alphaPos;
					f2 = Element.resSn+Element.ek*(d-Element.sn);
					f1 = Element.fP+Element.ekunload*deltaD ; 
					f = min(f1,f2);														
                elseif (d>xDevPos) 
					Element.ek = Element.capSlope*Element.elstk;
					f2 = yDevPos+Element.ek*(d-xDevPos);
					f1 = Element.fP+Element.ekunload*deltaD  ;
					f = min(f1,f2); 
                    if (Element.ek~=Element.ekunload) 
                         [f,Element]=envHitsZero(f,Element);
                    end					
                else
                end                
			elseif (Element.cpPos > Element.dyPos) 				
				[xDevPos1,yDevPos1]=interPoint(Element.sn,Element.resSn,Element.ekhardPos,Element.cpPos,Element.fCapPos,Element.capSlope*Element.elstk);                
				if(d<=Element.sn && d<=xDevPos1) % Unloading in positive loading direction			                
                    Element.ek = Element.ekunload;
					f  = Element.fP+Element.ek*deltaD;					
				elseif ((d>Element.sn)&&(d<=Element.cpPos) && d<=xDevPos1) 				
                    Element.ek = Element.elstk*Element.alphaPos;
					f2 = Element.resSn+Element.ek*(d-Element.sn);
					f1 = Element.fP+Element.ekunload*deltaD;
					f = min(f1,f2); 
                    if (abs(f-f1)<1.0e-10) 
                        Element.ek=Element.ekunload;
                    end					
				elseif((d>Element.cpPos)&&(d<Element.fracDispPos)) 
					Element.ek = Element.capSlope*Element.elstk;
					f2 = Element.fCapPos+Element.ek*(d-Element.cpPos);
					f1 = Element.fP+Element.ekunload*deltaD;  
					f = min(f1,f2); 
                    if (abs(f-f1)<1.0e-10) 
                        Element.ek=Element.ekunload;
                    end
                    if (Element.ek~=Element.ekunload) 
                        [f,Element]=envHitsZero(f,Element);
                    end					
                %added by Dimitrios to simulate ductile fracture
	            elseif(d>=Element.fracDispPos)
	                f=1.0e-10;
	                Element.ek = -1.0e-10;	                
	                Element.flagstopdeg = 1;
                end
                % Added by Dimitrios to stay on the residual path when  dresidual <d< dfracture               
                if(d < Element.fracDispPos &&...
                       d > Element.cpPos+(Element.Resfac*Element.fyieldPos-Element.fCapPos)/(Element.capSlope*Element.elstk))
                    f = Element.Resfac*Element.fyieldPos;
                    Element.ek = -1.0e-10;                    
                end
                if( Element.flagControlResponse == 1) % Response has to be zero since post capping regions hit zero
                    f = 1.0e-10;                    
                end                
            end

        % IF Element.LP IS EQUAL TO 1
        elseif (Element.LP==1) 
            if(d<=Element.sn) 
                Element.ek = Element.ekunload;
				f  = Element.fP+Element.ek*deltaD;				
            elseif (((d>Element.sn)&&(Element.sn==Element.snEnv)&&(d<=Element.snHor))||((Element.iNoFneg==1)&&(d>Element.sn)&&(d<Element.snHor))) 
				Element.ek = Element.elstk*Element.alphaPos;
				f2 = Element.resSn+Element.ek*(d-Element.sn);
				f1 = Element.fP+Element.ekunload*deltaD;  
				f = min(f1,f2); 
                if (abs(f-f1)<1.0e-10) 
                    Element.ek=Element.ekunload;
                end				
			else
				Element.ek = 0;
				f1 = Element.fP+Element.ekunload*deltaD;  
				f2 = Element.fLimPos;
				f = min(f1,f2); 
                if (abs(f-f1)<1.0e-10) 
                    Element.ek=Element.ekunload;
                end				
            end
        else
        end
    % Elastic
    else
        if (d>0.0) 
            [d,f,Element.ek]=envelPosCap2(Element.fyPos,Element.alphaPos,Element.capSlope,Element.cpPos,d,0.0,Element.ek,Element.elstk,Element.fyieldPos,Element.Resfac);
        else
            [d,f,Element.ek]=envelNegCap2(Element.fyNeg,Element.alphaNeg,Element.capSlopeNeg,Element.cpNeg,d,0.0,Element.ek,Element.elstk,Element.fyieldNeg,Element.ResfacNeg);
        end                       
    end
	
% IF   D E L T A < 0 - - - - - - - - - - - - - - - - - - - - - - - -  
else
    if (Element.iNoFneg==1) 
        [Element.dNewLoadNeg]=interPoint(Element.dyPos,Element.fyPos,Element.ekhardPos,0.0,0.0,0.0);
    end
    % If there is a corner changing delta from positive to negative ---      
    if (Element.kon==1) 
        Element.kon = 2;
		Element.dlstPos = Element.dP;
		Element.flstPos = Element.fP;

        if (Element.dlstPos>=Element.dmax) 
			Element.fLimPos = Element.flstPos;
			Element.dLimPos = Element.dlstPos;
        end

        if (Element.resSp<0) 
			Element.RSE = 0.5*Element.fP*Element.fP/Element.ekunload;
		else
			Element.RSE=0.5*(Element.fP+Element.resSn)*(Element.dP-Element.sp)+0.5*Element.resSp*Element.resSp/(Element.elstk*Element.alphaNeg);
        end	    
        
        % Update the Unloading Stiffness deterioration
        if((Element.ecapk~=0.0)&&(d>Element.sn)&&(Element.flagstopdeg~=1)&&(Element.flagdeg~=1)) 
            betak = ((Element.Enrgc-Element.RSE)/(Element.Enrgtk-(Element.Enrgtot-Element.RSE)))^Element.ck;
            if(((Element.Enrgtot-Element.RSE)>=Element.Enrgtk)||(betak>=1.0)) 
                betak = 1.0;
            end
            % If Post caping slopes have not been flat due to residual update stiffness deterioration
            if(Element.flagdeg ~=1 || Element.flagstopdeg~=1)	        
                Element.ekunload = Element.ekexcurs*(1-betak);
            else % Keep the same unloading stiffness
                Element.ekunload = Element.ekexcurs;
            end            
            if(Element.ekunload<=0.1*Element.elstk) 
                Element.ekunload = 0.1*Element.elstk;
            end
        end

        % Element.sp CALCULATION----------------------------------------------------
        if((Element.dmin<Element.dyieldNeg)||(Element.dmax>Element.dyieldPos)) 
            if((Element.dP>Element.sn)||(((Element.dP<Element.sn)&&(Element.ekP==Element.ekunload))))
                [Element]=spCalc(Element);
                if((abs(Element.dmin-Element.dyieldNeg)>=1.0e-10)&&(abs(Element.sp-Element.dyieldNeg)<=1.0e-10)) 
                    Element.sp=Element.dyieldNeg-1.0e-9;
                end
            end
        end

        if (Element.sp<Element.dmin) 
            Element.dmin = Element.sp;
        end

        if ((Element.iNoFpos==1)&&(Element.dP>=Element.dNewLoadPos)) 
            Element.sp = Element.dNewLoadPos + 1.0e-6;
            Element.resSp = 0;
        end
    end
    
    %	UNLOADING
    %  Push envelope
    if ((Element.iNoFneg==1)&&(Element.iNoFpos==1)&&(d>Element.fracDispNeg)) 
        f = Element.ResfacNeg*Element.fyieldNeg;
		Element.ek = 0;		
    elseif ((Element.iNoFneg==1)&&(Element.iNoFpos==1)&&(d<=Element.fracDispNeg)...
             ||Element.flagControlResponse==1)
        f = -1.0e-10;
	    Element.ek =0.0;	    
		Element.flagstopdeg = 1;
	elseif ((Element.iNoFneg==1)&&(d<Element.sp)&&(d>Element.fracDispNeg)) 
		f = Element.ResfacNeg*Element.fyieldNeg;
		Element.ek = 0.0;	
	elseif ((Element.iNoFneg==1)&&(d<Element.sp)&&(d<=Element.fracDispNeg)) 
		f = -1.0e-10;
		Element.ek = 0.0;		
		Element.flagstopdeg = 1;
	elseif (d<=Element.fracDispNeg || Element.flagControlResponse ==1)
		f = -1.0e-10;
		Element.ek = 0.0;		
	    Element.flagstopdeg = 1;
        Element.flagControlResponse = 1;               % To control response after I exceed fracture
	elseif ((Element.iNoFneg==1)&&(d>=Element.sp)) 
		Element.ek = Element.ekunload;
		f  = Element.fP+Element.ek*deltaD;		
	elseif ((Element.iNoFpos==1)&&(d>Element.dNewLoadPos)) 
		f = 0;
		Element.ek = 0;		
	elseif (d<Element.dmin) 
        [d,f,Element.ek]=envelNegCap2(Element.fyNeg,Element.alphaNeg,Element.capSlopeNeg,Element.cpNeg,d,0.0,Element.ek,Element.elstk,Element.fyieldNeg,Element.ResfacNeg);
        Element.dmin = d;		                      

        % COMPUTE MINIMUM POSSIBLE DISPLACEMENT   
        [dBoundNeg,Element]=boundNeg(Element);
        if((d<dBoundNeg)||(Element.ek==1.0e-7)) 
            Element.iNoFneg = 1;
        end
    elseif (abs(Element.sp)>1.0e-10)
        if (Element.LN==0) 
            if (Element.cpNeg>=Element.dyNeg) 
				[xDevNeg1,yDevNeg1]=interPoint(Element.sp,Element.resSp,Element.elstk*Element.alphaNeg,Element.cpNeg,Element.fCapNeg,Element.capSlopeNeg*Element.elstk);
				[xDevNeg2,yDevNeg2]=interPoint(Element.sp,Element.resSp,Element.ekunload,Element.cpNeg,Element.fCapNeg,Element.capSlopeNeg*Element.elstk);
                if (xDevNeg1<xDevNeg2) 
					xDevNeg = xDevNeg1;
					yDevNeg = yDevNeg1;
				else
					xDevNeg = xDevNeg2;
					yDevNeg = yDevNeg2;
                end	
                if ((d>=Element.sp)&&(d>=xDevNeg)) 
					Element.ek = Element.ekunload;
					f  = Element.fP+Element.ek*deltaD;					
                elseif ((d<Element.sp)&&(d>=xDevNeg)) 
					Element.ek = Element.elstk*Element.alphaNeg;
					f2 = Element.resSp+Element.ek*(d-Element.sp);
					f1 = Element.fP+Element.ekunload*deltaD;
					f = max(f1,f2);
                    if(Element.fyNeg >= Element.ResfacNeg*Element.fyieldNeg)                   
                        f  = Element.ResfacNeg*Element.fyieldNeg; 
                    end
                    if (abs(f-f1)<1.d-10) 
                        Element.ek=Element.ekunload;
                    end					
                elseif (d<xDevNeg) 
					Element.ek = Element.capSlopeNeg*Element.elstk;
					f2 = yDevNeg+Element.ek*(d-xDevNeg);
					f1 = Element.fP+Element.ekunload*deltaD; 
					f = max(f1,f2);
                    if (abs(f-f1)<1.0e-10) 
                        Element.ek=Element.ekunload;
                    end
                    if (Element.ek~=Element.ekunload) 
                        [f,Element]=envHitsZero(f,Element);
                    end					
                else
                end
			elseif (Element.cpNeg<Element.dyNeg) 			
				[xDevNeg1,yDevNeg1]=interPoint(Element.sp,Element.resSp,Element.elstk*Element.alphaNeg,Element.cpNeg,Element.fCapNeg,Element.	capSlopeNeg*Element.elstk);
				if ( d>=Element.sp && d>=xDevNeg1) 			                
                    Element.ek = Element.ekunload;
					f = Element.fP+Element.ek*deltaD;					
				elseif((d<Element.sp)&&(d>=Element.cpNeg)&&(d>=xDevNeg1)) 
					Element.ek = Element.elstk*Element.alphaNeg;
					f2 = Element.resSp+Element.ek*(d-Element.sp);
					f1 = Element.fP+Element.ekunload*deltaD;  
					f=max(f1,f2);
                    if (abs(f-f1)<1.d-10) 
                        Element.ek=Element.ekunload;
                    end					
                elseif((d<Element.cpNeg)&&(d>Element.fracDispNeg)) 
					Element.ek = Element.capSlopeNeg*Element.elstk;
					f2 = Element.fCapNeg+Element.ek*(d-Element.cpNeg);
 					f1 = Element.fP+Element.ekunload*deltaD;  
					f = max(f1,f2);
                    if (abs(f-f1)<1.d-10) 
                        Element.ek=Element.ekunload;
                    end
                    if (Element.ek~=Element.ekunload)
                        [f,Element]=envHitsZero(f,Element);
                    end					
                % Added by Dimitris to model ductile tearing. Once theta_u(fracDispNeg) is exceeded Strength drops to zero and stays
	            elseif(d<=Element.fracDispNeg ||  Element.flagControlResponse == 1)
	 			    Element.ek = -1.0e-10;
	                f = 1.0e-10;	                
					Element.flagstopdeg = 1;
                    Element.flagControlResponse = 1;               % To dictate the response after passing fracture             
                end
                % Added by Dimitrios to stay on the residual path when  dresidual <d- < dfracture               
                if(d > Element.fracDispNeg &&...
                        d < Element.cpNeg+(Element.ResfacNeg*Element.fyieldNeg-Element.fCapNeg)/(Element.capSlopeNeg*Element.elstk))
                    f = Element.ResfacNeg*Element.fyieldNeg;
                    Element.ek = -1.0e-10;                    
                end
                if( Element.flagControlResponse == 1) % Response has to be zero since post capping regions hit zero
                    f = 1.0e-10;                    
                end
            end
		elseif (Element.LN==1)
            if(d>=Element.sp)
				Element.ek = Element.ekunload;
				f  = Element.fP+Element.ek*deltaD;				
			elseif (((d<Element.sp)&&(Element.sp==Element.spEnv)&&(d>Element.spHor))||((Element.iNoFpos==1)&&(d<Element.sp)&&(d>Element.spHor))) 
				Element.ek = Element.elstk*Element.alphaNeg;
				f2 = Element.resSp+Element.ek*(d-Element.sp);
				f1 = Element.fP+Element.ekunload*deltaD;
				f = max(f1,f2);
                if (abs(f-f1)<1.0e-10)
                    Element.ek=Element.ekunload;
                end				
            else
				Element.ek = 0;
				f1 = Element.fP+Element.ekunload*deltaD ; 
				f2 = Element.fLimNeg;
				f = max(f1,f2);
                if (abs(f-f1)<1.d-10) 
                    Element.ek=Element.ekunload;
                end				
            end
        else
        end
    else
        if (d>0.0) 
            [d,f,Element.ek]=envelPosCap2(Element.fyPos,Element.alphaPos,Element.capSlope,Element.cpPos,d,0.0,Element.ek,Element.elstk,Element.fyieldPos,Element.Resfac);
        else
            [d,f,Element.ek]=envelNegCap2(Element.fyNeg,Element.alphaNeg,Element.capSlopeNeg,Element.cpNeg,d,0.0,Element.ek,Element.elstk,Element.fyieldNeg,Element.ResfacNeg);
        end                
    end
end
%%	ENDS BIG LOOP ****************************************************
% 	       
Enrgi = 0.5*(f+Element.fP)*deltaD;
Element.Enrgc = Element.Enrgc + Enrgi;
Element.Enrgtot = Element.Enrgtot + Enrgi;
Element.RSE = 0.5*f*f/Element.ekunload;

%	Flag to deteriorate parameters on the opposite side of the loop --	
if((f*Element.fP<0.0)&&(Element.interup==0)) 
    if(((Element.fP>0.0)&&(Element.dmax>Element.dyieldPos))||((Element.fP<0.0)&&(Element.dmin<Element.dyieldNeg))) 
        Element.flagdeg = 1;		
        Element.interup=1;
    end
end	


%	Element.energy CALCULATIONS ---------------------------------------------
if((Element.flagstopdeg==0)&&(Element.flagdeg==1))             	
    if((Element.Enrgtot>=Element.Enrgts)&&(Element.Enrgts~=0.0))
        betas = 1.0;
    elseif((Element.Enrgtot>=Element.Enrgtd)&&(Element.Enrgtd~=0.0)) 
        betad = 1.0;
    else
         if(Element.ecaps~=0.0) 
            betas = (Element.Enrgc/(Element.Enrgts-Element.Enrgtot))^Element.cs;
         end
         if(Element.ecapd~=0.0) 
            betad = (Element.Enrgc/(Element.Enrgtd-Element.Enrgtot))^Element.cd;
         end

         if(abs(betas)>=1.0) 
            betas = 1.0;
         end
         if(abs(betad)>=1.0) 
            betad = 1.0;
         end
    end 
    %  Initialize Element.energy of the cycle and Kstif for next loop -------
    Element.Enrgc = 0.0;
    Element.ekexcurs = Element.ekunload;   

    % Deteriorate parameters for the next half cycle
    if(deltaD<0.0)
        if(Element.flagstopdeg == 0)
			Element.fyNeg = Element.fyNeg*(1-betas*Element.DNeg);
			Element.alphaNeg=Element.alphaNeg*(1-betas*Element.DNeg);
			Element.fCapRefNeg=Element.fCapRefNeg*(1-betad*Element.DNeg);		
		end
		% When we reach post capping slope goes to zero due to residual
		if(Element.fyNeg>=Element.ResfacNeg*Element.fyieldNeg) % If strength drops below residual
			Element.fyNeg = Element.ResfacNeg*Element.fyieldNeg;
			Element.alphaNeg = 10^(-4);
			Element.fCapRefNeg = Element.fyNeg;
			Element.capSlopeNeg = -10^(-6);
			Element.flagstopdeg = 1;
		else % Keep updating the post capping slope
			Element.capSlopeNeg = Element.capSlopeOrigNeg*(1-abs((Element.ResfacNeg*Element.fyieldNeg)/Element.fyNeg));
			if(Element.capSlopeNeg >=0)
				Element.capSlopeNeg = -10^(-6);
			end
		end

        Element.dyNeg = Element.fyNeg/Element.elstk;
        Element.ekhardNeg=Element.alphaNeg*Element.elstk;

        Element.dCap1Neg=Element.fCapRefNeg/(Element.elstk-Element.capSlopeOrigNeg*Element.elstk);
        Element.dCap2Neg=(Element.fCapRefNeg+Element.ekhardNeg*Element.dyNeg-Element.fyNeg)/(Element.ekhardNeg-Element.capSlopeOrigNeg*Element.elstk);
        Element.cpNeg=min(Element.dCap1Neg,Element.dCap2Neg);

        Element.fCapNeg = Element.fCapRefNeg + Element.capSlopeOrigNeg*Element.elstk*Element.cpNeg;

        [Element.dLimNeg,Element.fLimNeg,Element.ekt]=envelNegCap2(Element.fyNeg,Element.alphaNeg,Element.capSlopeNeg,Element.cpNeg,Element.dLimNeg,Element.fLimNeg,Element.ekt,Element.elstk,Element.fyieldNeg,Element.ResfacNeg);
        [Element]=spCalc(Element);
        % 	In case the degradation point moves from the negative to pos. side
        if(Element.resSp>f)         
            d = Element.sp;
            f = Element.resSp;
            Element.ek = Element.ekhardNeg;   
        end
    else
        if(Element.flagstopdeg == 0)
			Element.fyPos = Element.fyPos*(1-betas*Element.DPlus);	
			Element.alphaPos=Element.alphaPos*(1-betas*Element.DPlus);	
			Element.fCapRefPos=Element.fCapRefPos*(1-betad*Element.DPlus);		
		end
                
        %If post capping slope goes to zero due to residual:
        if(Element.fyPos <= Element.Resfac*Element.fyieldPos) % If yield Strength Pos drops below residual
			Element.fyPos = Element.Resfac*Element.fyieldPos;
			Element.alphaPos = 10^(-4);
			Element.fCapRefPos = Element.fyPos;
			Element.capSlope = -10^(-6);
			Element.flagstopdeg = 1;              
        else % keep updating
		    Element.capSlope = Element.capSlopeOrig*(1-abs((Element.Resfac*Element.fyieldPos)/Element.fyPos));
            if(Element.capSlope >=0)
                Element.capSlope = -10^(-6);
            end
        end
        Element.dyPos = Element.fyPos/Element.elstk;
        Element.ekhardPos=Element.alphaPos*Element.elstk;

        Element.dCap1Pos=Element.fCapRefPos/(Element.elstk-Element.capSlopeOrig*Element.elstk);
        Element.dCap2Pos=(Element.fCapRefPos+Element.ekhardPos*Element.dyPos-Element.fyPos)/(Element.ekhardPos-Element.capSlopeOrig*Element.elstk);
        Element.cpPos=max(Element.dCap1Pos,Element.dCap2Pos);
        Element.fCapPos = Element.fCapRefPos + Element.capSlopeOrig*Element.elstk*Element.cpPos;


        [Element.dLimPos,Element.fLimPos,Element.ekt]=envelPosCap2(Element.fyPos,Element.alphaPos,Element.capSlope,Element.cpPos,Element.dLimPos,Element.fLimPos,Element.ekt,Element.elstk,Element.fyieldPos,Element.Resfac);
        [Element]=snCalc(Element);
        %In case the degradation point moves from the pos. to neg. side
        if(Element.resSn<f) 
            d = Element.sn;
            f = Element.resSn;
            Element.ek = Element.ekhardPos;
        end	
    end
end

% Check the horizontal limit in case that dBound is reached after first neg slope
if ((d<0)&&(abs(Element.ek)<=1.0d-7)) 
    Element.LN = 1;
end
if ((d>0)&&(abs(Element.ek)<=1.0e-7)) 
    Element.LP = 1;
end

% Calculation of static resisting force vector for the element -----
RestoringForce=zeros(6,1);
RestoringForce(3,1)=f;
RestoringForce(6,1)=-f;

if( Structure.RigidLinkNo )
    RestoringForce = TR' * RestoringForce;  
end
	
% Updating parameters for next cycle ---------------------------------
Element.ekP = Element.ek;
Element.fP = f;
%Element.dP = d;

%priority of logical operators
if (Element.interup==1&&(Element.ek==Element.elstk*Element.alphaPos)...
                        ||(Element.ek==Element.elstk*Element.alphaNeg)...
                        ||(Element.ek==Element.capSlope*Element.elstk)...
                        ||(Element.ek==Element.capSlopeNeg*Element.elstk)) 
    Element.interup = 0;
end
	




%% BEGIN FUNCTIONS %%
function [Element]=snCalc(Element)
% c	Each time that the hysteretic loop changes from negative to positive
% c	delta displacement, this subroutine calculates the point where ends
% c	Element.ekunload and starts	whether the positive hardening curve or the 
% c	positive cap curve. In cases where "Element.cpPos<Element.fyPos" this point could
% c	not be reached, because	the unloading stiffness could intersect the
% c	positive cap curve. However, this change is reflected in the main 
% c	program.
% c	
% c	Output Variables: Element.sn,Element.resSn,snHard,resSnHard
% c	Input Variables:  Element.dP,Element.fP,Element.ekunload,Element.alphaPos,Element.dyPos,Element.fyPos,Element.cpPos,Element.fCapPos
% c	     			  capStiff,Element.fCapRefPos

Resid = Element.Resfac*Element.fyPos;
dresid = Element.cpPos+(Resid-Element.fCapPos)/(Element.capSlope*Element.elstk);
ekresid = 1.0e-10;
Element.dyPos = Element.fyPos/Element.elstk;

if (Element.dyPos<Element.cpPos) 
    [snHard,resSnHard]=interPoint(Element.dyPos,Element.fyPos,Element.elstk*Element.alphaPos,Element.dP,Element.fP,Element.ekunload);
else
    [snHard,resSnHard]=interPoint(Element.cpPos,Element.fCapPos,Element.elstk*Element.alphaPos,Element.dP,Element.fP,Element.ekunload);
end

[snCap,resSnCap]=interPoint(0.0,Element.fCapRefPos,Element.capSlope*Element.elstk,Element.dP,Element.fP,Element.ekunload);

Element.sn = min(snHard,snCap);
Element.resSn = min(resSnHard,resSnCap);
Element.snEnv = Element.sn;

if((Element.LP==1)&&(Element.fLimPos==0.0)) 
    [snLim,resSnLim]=interPoint(Element.dLimPos,Element.fLimPos,0.d0,Element.dP,Element.fP,Element.ekunload);
    if (snLim<Element.sn) 
        Element.sn=snLim;
        Element.resSn=resSnLim;
    end

    [Element.snHor]=interPoint(Element.dLimPos,Element.fLimPos,0.d0,Element.dyPos,Element.fyPos,Element.elstk*Element.alphaPos);
end

if (Element.sn>dresid) 
    [snResid,resSnResid]=interPoint(dresid,Resid,ekresid,Element.dP,Element.fP,Element.ekunload);
    Element.sn = snResid;
    Element.resSn = resSnResid;
end
	
% c	*******************************************************************
% 	subroutine spCalc(Element.sp,Element.resSp,Element.spEnv,Element.resSpEnv,Element.dP,Element.fP,Element.ekunload,
%      &				Element.alphaNeg,Element.dyNeg,Element.fyNeg,Element.cpNeg,Element.fCapNeg,Element.capSlope,
%      &		Element.fCapRefNeg,Element.LN,Element.dLimNeg,Element.fLimNeg,Element.spHor,Element.resSpHor,Element.elstk,Element.Resfac)
% 
% c	Idem as snCalc in the opposite direction.
% c
% c	Output Variables: Element.sp,Element.resSp
% c	Input Variables:  Element.dP,Element.fP,Element.ekunload,Element.alphaNeg,Element.dyNeg,Element.fyNeg,Element.cpNeg,Element.fCapNeg
% c	     			  capStiff,Element.fCapRefNeg
% c	*******************************************************************
% 
% c	LABELLED COMMONS
% c      include 'infel13.h'
% 
% 	integer Element.LN
% 
% 	real*8 Element.sp,Element.resSp,Element.dP,Element.fP,Element.ekunload,Element.elstk,Element.alphaNeg,Element.dyNeg,
%      &	   Element.fyNeg,Element.Resfac,Resid,resSpCap,
%      &	   Element.fCapNeg,dresid,Element.capSlope,Element.cpNeg,spCap,Element.fCapRefNeg,ekresid,
%      &	   spResid,resSpResid,Element.spEnv,Element.resSpEnv,Element.spHor,Element.resSpHor,spLim,
%      &	   resSpLim,Element.dLimNeg,Element.fLimNeg,spHard,resSpHard
% 
function [Element]=spCalc(Element)
        
Resid = Element.ResfacNeg*Element.fyNeg;
Element.dyNeg = Element.fyNeg/Element.elstk;
dresid = Element.cpNeg+(Resid-Element.fCapNeg)/(Element.capSlopeNeg*Element.elstk);
ekresid = 1.0e-10;

if (Element.dyNeg>Element.cpNeg) 
    [spHard,resSpHard]=interPoint(Element.dyNeg,Element.fyNeg,Element.elstk*Element.alphaNeg,Element.dP,Element.fP,Element.ekunload);
else
    [spHard,resSpHard]=interPoint(Element.cpNeg,Element.fCapNeg,Element.elstk*Element.alphaNeg,Element.dP,Element.fP,Element.ekunload);
end

[spCap,resSpCap]=interPoint(0.0,Element.fCapRefNeg,Element.capSlopeNeg*Element.elstk,Element.dP,Element.fP,Element.ekunload);
Element.sp = max(spHard,spCap);
Element.resSp = max(resSpHard,resSpCap);
Element.spEnv = Element.sp;

if((Element.LN==1)&&(Element.fLimNeg==0.0))
    [spLim,resSpLim]=interPoint(Element.dLimNeg,Element.fLimNeg,0.d0,Element.dP,Element.fP,Element.ekunload);
    if (spLim>Element.sp) 
        Element.sp=spLim;
        Element.resSp=resSpLim;
    end

    [Element.spHor]=interPoint(Element.dLimNeg,Element.fLimNeg,0.0,Element.dyNeg,Element.fyNeg,Element.elstk*Element.alphaNeg);
end

if (Element.sp<dresid) 
    [spResid,resSpResid]=interPoint(dresid,Resid,ekresid,Element.dP,Element.fP,Element.ekunload);
    Element.sp = spResid;
    Element.resSp = resSpResid;
end	

% 
% 
% c	*******************************************************************
% 	subroutine boundPos (dBoundPos,Element.fCapRefPos,Element.capSlope,Element.fyNeg,Element.dyNeg,
%      &				Element.alphaNeg,Element.fCapPos,Element.cpPos,Element.elstk)	
% 	
% c	Each time that maximum displacement is exceeded,this subroutine
% c	calculates the maximum displacement that can occur on the positive
% c	side. The max. displacement is based on the intersection of the
% c	current path with the extension of the negative hardening curve
% c	(i.e. opposite side)
% c
% c	Output Variables: dBoundPos
% c	Input Variables: Element.fCapRefPos,capStiff,Element.fyNeg,Element.dyNeg,Element.alphaNeg,Element.cpPos,
% c				     Element.fCapPos	
% c	*******************************************************************
% 
function [dBoundPos,Element]=boundPos(Element)
Resid=Element.Resfac*Element.fyieldPos;    

Element.dyNeg = Element.fyNeg/Element.elstk;
dresid = Element.cpPos+(Resid-Element.fCapPos)/(Element.capSlope*Element.elstk);
ekresid = 1.0e-10;

[d1]=interPoint(Element.dyNeg,Element.fyNeg,Element.elstk*Element.alphaNeg,0.0,Element.fCapRefPos,Element.capSlope*Element.elstk);

[d2]=interPoint(Element.dyNeg,Element.fyNeg,Element.elstk*Element.alphaNeg,dresid,Resid,ekresid);
dBoundPos = max(d1,d2);

	
%%
% 
% c	*******************************************************************
% 	subroutine boundNeg (dBoundNeg,Element.fCapRefNeg,Element.capSlope,Element.fyPos,Element.dyPos,
%      &				Element.alphaPos,Element.fCapNeg,Element.cpNeg,Element.elstk)	
% c
% c	Idem as boundPos but for the negative case.
% c
% c	Output Variables: dBoundNeg
% c	Input Variables: Element.fCapRefNeg, capStiff,Element.fyPos,Element.dyPos,Element.alphaPos,Element.cpNeg,
% c					 Element.fCapNeg	
% c	*******************************************************************

function [dBoundNeg,Element]=boundNeg(Element)

Resid = Element.ResfacNeg*Element.fyieldNeg;
    
Element.dyPos = Element.fyPos/Element.elstk;
dresid = Element.cpNeg+(Resid-Element.fCapNeg)/(Element.capSlopeNeg*Element.elstk);
ekresid = 1.0e-10;

[d1]=interPoint(Element.dyPos,Element.fyPos,Element.elstk*Element.alphaPos,0.d0,Element.fCapRefNeg,Element.capSlopeNeg*Element.elstk);
[d2]=interPoint(Element.dyPos,Element.fyPos,Element.elstk*Element.alphaPos,dresid,Resid,ekresid);
dBoundNeg = min(d1,d2);

%%
% c	*******************************************************************
%  	subroutine envHitsZero(f,Element.fP,Element.ek)
% c	
% c	When the horizontal axis is crossed with the loading path on the 
% c	Element.capSlope, the force is taken as zero and a flag is displayed
% c Comments Added By DLignos
% c For positive path crossing zero: Flag: iNoFpos = 1 --> When crosses Zero Force
% c For negative path crossing zero: Flag: iNoFneg = 1 --> When crosses Zero Force
% c If The element hits zero from negative or positive loading direction
%   must stay at zero Strength:
% c Element.flagControlResponse = 1;
% c	*******************************************************************
% 
% 	integer Element.iNoFneg,Element.iNoFpos,NoCollapse,maybeStop
% 	real*8 f,Element.fP,Element.ek
% 	common/colapso/Element.iNoFneg,Element.iNoFpos,NoCollapse,maybeStop
function [f,Element]=envHitsZero(f,Element)
if (Element.fP>0) 
    if ((f*Element.fP)<0) 
        f = 0;
    Element.ek = 0;
        Element.iNoFpos = 1;
        Element.flagControlResponse = 1;
    end
elseif (Element.fP<0) 
    if ((f*Element.fP)<0) 
        f = 0;
        Element.ek = 0;
        Element.iNoFneg = 1;
        Element.flagControlResponse = 1;
    end
else
end

% 	subroutine interPoint(xInt,yInt,x1,y1,m1,x2,y2,m2)
function [xInt,yInt]=interPoint(x1,y1,m1,x2,y2,m2)
xInt = (-m2*x2+y2+m1*x1-y1) / (m1-m2);
yInt = m1*xInt-m1*x1+y1;


%	subroutine envelPosCap2(fy,alphaPos,alphaCap,cpDsp,d,f,ek,elstk,fyieldPos,Resfac)
function [d,f,ek]=envelPosCap2(fy,alphaPos,alphaCap,cpDsp,d,f,ek,elstk,fyieldPos,Resfac)
fracDispPos=0.20;
dy = fy/elstk;

if(dy<=cpDsp) 
    Res = Resfac*fyieldPos;
    rcap = fy+alphaPos*elstk*(cpDsp-dy);
    dres = cpDsp+(Res-rcap)/(alphaCap*elstk);

    if (d<0.0)
        f = 0.0;
        ek = 0.0;
    elseif (d<=dy) 
        ek = elstk;
        f = ek*d;
    elseif(d<=cpDsp) 
        ek = elstk*alphaPos;
        f = fy+ek*(d-dy);
    elseif(d<=dres) 
        ek = alphaCap*elstk;
        f = rcap+ek*(d-cpDsp);
    else
        ek = 1.0e-7;
        f = Res+d*ek;
    end
    % added by Dimitrios to account for fracture	
    if(d>=fracDispPos)
        ek = 1.0e-7;
        f = 1.0e-10;
        d=fracDispPos;        
    end
elseif(dy>cpDsp) 
    rcap = elstk*cpDsp;
    Res = Resfac*rcap;
    dres = cpDsp+(Res-rcap)/(alphaCap*elstk);

    if (d<0.0) 
        f = 0.0;
        ek = 0.0;
    elseif(d<=cpDsp)
        ek = elstk;
        f = ek*d;
    elseif(d<=dres)
        ek = alphaCap*elstk;
        f = rcap+ek*(d-cpDsp);
    else
        ek = 1.0e-7;
        f = Res+d*ek;
    end
    % added by Dimitrios to account for fracture	
    if(d>=fracDispPos)
        ek = 1.0e-7;
        f = 1.0e-10;
        d=fracDispPos;        
    end
else
end

% 	subroutine envelNegCap2(fy,alphaNeg,alphaCap,cpDsp,d,f,ek,elstk,fyieldNeg,Resfac)
function [d,f,ek]=envelNegCap2(fy,alphaNeg,alphaCap,cpDsp,d,f,ek,elstk,fyieldNeg,Resfac)

fracDispNeg = -0.20;
dy = fy/elstk;

if(dy>=cpDsp) 
    Res = Resfac*fyieldNeg;
    rcap = fy+alphaNeg*elstk*(cpDsp-dy);
    dres = cpDsp+(Res-rcap)/(alphaCap*elstk);

    if (d>0.0) 
        f = 0.0;
        ek = 0.0;
    elseif (d>=dy) 
        ek = elstk;
        f = ek*d;
    elseif (d>=cpDsp)
        ek = elstk*alphaNeg;
        f = fy+ek*(d-dy);
    elseif (d>=dres)
        ek = elstk*alphaCap;
        f = rcap+ek*(d-cpDsp);
    else
        ek = 1.0e-7;
        f = Res+ek*d;
    end
    % added by Dimitrios to account for fracture	
    if(d<=fracDispNeg)
        ek = 1.0e-7;
        f = 1.0e-10;
        d=fracDispNeg;        
    end
	
elseif(dy<cpDsp)     
    rcap = elstk*cpDsp;
    Res = Resfac*rcap;
    dres = cpDsp+(Res-rcap)/(alphaCap*elstk);

    if (d>0.0) 
        f = 0.0;
        ek = 0.0;
    elseif (d>=cpDsp)
        ek = elstk;
        f = ek*d;
    elseif (d>=dres)
        ek = elstk*alphaCap;
        f = rcap+ek*(d-cpDsp);
    else
        ek = 1.0e-7;
        f = Res+ek*d;
    end
    % added by Dimitrios to account for fracture	
    if(d<=fracDispNeg)
        ek = 1.0e-7;
        f = 1.0e-10;
        d=fracDispNeg;        
    end
else
end


