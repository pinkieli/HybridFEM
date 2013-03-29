%% Create Element Type 4
%Written by: Theodore L. Karavasilis
%This file constructs a SNAP (Stanford U) bilinear
%zero-length element with stiffenss and strength deterioration capabilities
%needs 29 input parameters  
function vars = CreateElement_Type5(ID, Type, Node1, Node2, DampStiffFac, DampMassFac, Ke, As, AsNeg, My_pos, My_neg, LamdaS, LamdaK, LamdaA, LamdaD, Cs, Ck, Ca, Cd,...
                              Thetap_pos, Thetap_neg, Thetapc_pos, Thetapc_neg, K, KNeg, Thetau_pos, Thetau_neg, DPlus, DNeg)
% Set element nodes
vars.ID = ID;
vars.Type=Type;           
vars.Nodes = [Node1 Node2];
vars.DampStiffFac=DampStiffFac;
vars.DampMassFac=DampMassFac;

% Set handles to element matrices forming function and restoring force
vars.FormElementMatrices = @FormElementMatrices_Type5;
vars.GetRestoringForce = @GetRestoringForce_Type5;

% keep track of element total displacement 
vars.Uprev=0.0;  

% initialize infel to zero
vars.kon=0; 
vars.dNewLoadPos=0.0; 
vars.dNewLoadNeg=0.0; 
vars.flagdeg=0; 
vars.flagstopdeg=0; 
vars.ekt=0; 
vars.interup=0.0;  
vars.iNoFneg=0.0; 
vars.iNoFpos=0.0; 
vars.LP=0.0; 
vars.LN=0.0;  
vars.dmax=0.0; 
vars.dmin=0.0; 
vars.Enrgtot=0.0; 
vars.fyPos=0.0; 
vars.fLimNeg=0.0; 
vars.fyNeg=0.0; 
vars.ekunload=0.0;  
vars.sp=0.0; 
vars.sn=0.0; 
vars.dP=0.0; 
vars.fP=0.0; 
vars.ek=0.0; 
vars.dLimPos=0.0; 
vars.dLimNeg=0.0; 
vars.cpPos=0.0; 
vars.cpNeg=0.0; 
vars.fLimPos=0.0; 
vars.dlstPos=0.0; 
vars.flstPos=0.0; 
vars.dlstNeg=0.0; 
vars.flstNeg=0.0; 
vars.ekexcurs=0.0; 
vars.RSE=0.0; 
vars.dCap1Pos=0.0; 
vars.dCap2Pos=0.0; 
vars.dCap1Neg=0.0; 
vars.dCap2Neg=0.0; 
vars.alphaNeg=0.0;  
vars.alphaPos=0.0; 
vars.ekhardNeg=0.0; 
vars.ekhardPos=0.0; 
vars.fCapRefPos=0.0; 
vars.fCapRefNeg=0.0; 
vars.Enrgts=0.0; 
vars.Enrgtk=0.0; 
vars.Enrgtd=0.0; 
vars.dyPos=0.0; 
vars.dyNeg=0.0; 
vars.dyieldPos=0.0; 
vars.dyieldNeg=0.0; 
vars.resSp=0.0; 
vars.resSn=0.0; 
vars.fCapPos=0.0; 
vars.fCapNeg=0.0; 
vars.snHor=0.0; 
vars.spHor=0.0; 
vars.snEnv=0.0; 
vars.spEnv=0.0; 
vars.capSlopeOrig=0.0; 
vars.Enrgc = 0.0;
vars.flagControlResponse = 0;    % Controls response after fracture

% set element properties from input
vars.elstk          = Ke;
vars.ekP            = Ke;
vars.fyieldPos      = My_pos;
vars.fyieldNeg      = My_neg;
vars.alpha          = As;
vars.alphaN         = AsNeg;
vars.ecaps          = LamdaS/(vars.fyieldPos/(vars.elstk));
vars.ecapk          = LamdaK/(vars.fyieldPos/(vars.elstk));
vars.ecapd		    = LamdaD/(vars.fyieldPos/(vars.elstk));
vars.cs		        = Cs;
vars.ck		        = Ck;
vars.cd		        = Cd;
vars.capDispPos     = Thetap_pos+vars.fyieldPos/vars.elstk;
vars.capDispNeg     = -Thetap_neg+vars.fyieldNeg/vars.elstk;
vars.Resfac         = K;
vars.ResfacNeg      = KNeg;                  
vars.capSlope	    = (-(vars.fyieldPos+vars.alpha*vars.elstk*vars.capDispPos))/(Thetapc_pos*vars.elstk);
vars.capSlopeNeg    = (-(abs(vars.fyieldNeg)+vars.alpha*vars.elstk*abs(vars.capDispNeg)))/(Thetapc_neg*vars.elstk); % added by DLignos
vars.fracDispPos    = Thetau_pos;
vars.fracDispNeg    = -Thetau_neg;
vars.DPlus          = DPlus;                   % For composite Action for Positive loading direction DPlus <= 1
vars.DNeg           = DNeg;                    % For composite Action for Negative loading direction DNeg <= 1
                                               % IF DPlus = DNeg = 1.0 then
                                               % simulation is without composite action