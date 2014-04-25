%% Create Material object
% data(1) = material ID
% data(2) = material type
%         = 1 : elastic 
%         = 2 : bilinear elasto plastic
%         = 3 : hysteretic material
% data(3:n) =  material properties 
function mat = CreateMaterial(data)

% Material structs need to be predefined with empty variables for code
% generation in Simulink
mat.E = 0;
mat.Et = 0;
mat.Fs = 0;
mat.Fs1 = 0;
mat.Fs2 = 0;
mat.Fy1 = 0;
mat.Fy2 = 0;
mat.ID = 0;
mat.K1 = 0;
mat.K2 = 0;
mat.K3 = 0;
mat.Kt = 0;
mat.Type = 0;
mat.YieldCode = 0;
mat.cd = 0;
mat.e = 0;
mat.ep = 0;
mat.es = 0;
mat.sb = 0;
mat.sc = 0;
mat.sy = 0;
mat.xprev = 0;
mat.c = 0;
mat.k1 = 0;
mat.k2 = 0;
mat.alpha = 0;
mat.uy = 0;
mat.a = 0;
mat.uref = 0;
mat.beta = 0;
mat.gamma1 = 0;
mat.n = 0;
mat.uref = 1e-24;
mat.Dmax = 0;   
mat.Dmin = 0; 
mat.zj = 0;
mat.Uprev = 0;
mat.k = 0;
mat.Dmax = 0;
mat.Dmin = 0;     
mat.Qj = 0;

mat.fpc = 0;
mat.epsc0 = 0;
mat.fpcu = 0;
mat.epscu = 0;
mat.CminStrain = 0;
mat.CendStrain = 0;
mat.Cstress = 0;
mat.Cstrain = 0;
mat.Ctangent = 0;
mat.CunloadSlope = 0;
mat.TminStrain = 0;
mat.TendStrain = 0;
mat.Tstrain = 0;
mat.Tstress = 0;
mat.Ttangent = 0;
mat.TunloadSlope = 0;

mat.Fy = 0;
mat.E0 = 0;
mat.b = 0;
mat.R0 = 0; 
mat.cR1 = 0;
mat.cR2 = 0;
mat.a1 = 0;   
mat.a2 = 0; 
mat.a3 = 0; 
mat.a4 = 0; 
mat.sigini = 0.; 
mat.Esh = 0.;
mat.epsy = 0.;
mat.konP = 0;
mat.kon = 0;
mat.eP = 0.;
mat.epsP = 0.0;
mat.sigP = 0.0;
mat.sig = 0.0;
mat.eps = 0.0;
mat.e = 0.;
mat.epsmax = 0.;
mat.epsmin = 0.0;
mat.epspl = 0.0;
mat.epss0 = 0.0;
mat.sigs0 = 0.0;
mat.epsr = 0.0;
mat.sigr = 0.0;
mat.epsmaxP = 0.;
mat.epsminP = 0.;
mat.epsplP = 0.0;
mat.epss0P = 0.0;
mat.sigs0P = 0.0;
mat.epssrP = 0.0;
mat.sigsrP = 0.0;


% Get material data line
material = str2num(data);

% Parse data line into ID and Type
mat.ID   = material(1);
mat.Type = material(2);   

% parse the material properties
switch (mat.Type)
    case 1  % elastic 
        mat.E     = material(3);			% modulus of elasticity
    case 2   % bilinear elasto plastic
        mat.e     = material(3);			% modulus of elasticity
        mat.sy    = material(4);			% yield stress
        mat.ep    = material(3)*material(5); % second slope
        mat.sc    = 0.;                     % current stress
        mat.sb    = 0.;                     % backstress  
        mat.cd    = 0;                      % yield code  
                                            %	=0, elastic
                                            %	=1, plastic
        mat.es    = 0.;                     % current strain
    case 3  % hysteretic material properties    
        mat.pinchX = material(3);
        mat.pinchY = material(4); 
        mat.damfc1 = material(5);
        mat.damfc2 = material(6); 
        mat.beta   = material(7);
        mat.mom1p  = material(8);
        mat.rot1p  = material(9);
        mat.mom2p  = material(10);
        mat.rot2p  = material(11);
        mat.mom3p  = material(12);
        mat.rot3p  = material(13);
        mat.mom1n  = material(14);
        mat.rot1n  = material(15);
        mat.mom2n  = material(16);
        mat.rot2n  = material(17);
        mat.mom3n  = material(18);
        mat.rot3n  = material(19);
        mat.energyA = 0.5*(mat.rot1p*mat.mom1p+ (mat.rot2p-mat.rot1p)*...
            (mat.mom2p+mat.mom1p)+ (mat.rot3p-mat.rot2p)*...
            (mat.mom3p+mat.mom2p)+ mat.rot1n*mat.mom1n +...
            (mat.rot2n-mat.rot1n)*(mat.mom2n+mat.mom1n)+...
            (mat.rot3n-mat.rot2n)*(mat.mom3n+mat.mom2n));

        mat.E1p = mat.mom1p/mat.rot1p;
        mat.E2p = (mat.mom2p-mat.mom1p)/(mat.rot2p-mat.rot1p);
        mat.E3p = (mat.mom3p-mat.mom2p)/(mat.rot3p-mat.rot2p);

        mat.E1n = mat.mom1n/mat.rot1n;
        mat.E2n = (mat.mom2n-mat.mom1n)/(mat.rot2n-mat.rot1n);
        mat.E3n = (mat.mom3n-mat.mom2n)/(mat.rot3n-mat.rot2n);

        mat.Eup = mat.E1p;
        if (mat.E2p> mat.Eup) mat.Eup = mat.E2p; end;
        if (mat.E3p> mat.Eup) mat.Eup = mat.E3p; end;

        mat.Eun = mat.E1n;
        if (mat.E2n > mat.Eun) mat.Eun = mat.E2n; end;
        if (mat.E3n > mat.Eun) mat.Eun = mat.E3n; end;

        mat.CrotMax  = 0.0;
        mat.CrotMin  = 0.0;
        mat.CrotPu   = 0.0;
        mat.CrotNu   = 0.0;
        mat.CenergyD = 0.0;
        mat.CloadIndicator = 0;
        mat.Cstress = 0.0;
        mat.Cstrain = 0.0;

        mat.TrotMax  = 0.0;
        mat.TrotMin  = 0.0;
        mat.TrotPu   = 0.0;
        mat.TrotNu   = 0.0;
        mat.TenergyD = 0.0;
        mat.TloadIndicator = mat.CloadIndicator;
        mat.Tstress = 0.0;
        mat.Tstrain = 0.0;
        mat.Ttangent = mat.E1p;
    case 4  % Bouc-Wen Model
        mat.c = material(3);      
        mat.k1 = material(4);    
        mat.k2= material(5);   
        mat.alpha= material(6);  
        mat.uy= material(7);    
        mat.a = material(8);
        mat.uref= 0.8*mat.uy;
        mat.beta= material(9);
        mat.gamma= material(10);
        mat.n = material(11);   
        mat.Dmax = 0.;   
        mat.Dmin = 0.; 
        mat.zj = 0.;
        mat.Uprev=0.;
        mat.k = mat.k1 + mat.k2;    
        mat.Dmax = 0.;
        mat.Dmin = 0.;        
    case 5  % Trilinear Model for panel zone shear deformation        
       % input material parameters 
       Vy = material(3);
       Vu = material(4);       
       Dy = material(5);
       Du = material(6);
       % equivalent spring stiffness
       Kinit = Vy/Dy;
       Ky    = (Vu-Vy)/(Du-Dy);       
       Ku    = 0.0;   
       % decompose spring stiffness 
	   mat.K1 = Kinit-Ky;
	   mat.K2 = Ky - Ku;
       mat.K3 = Ku;       
       mat.Kt = Kinit;
   	   mat.Fy1= mat.K1*Dy;
       mat.Fy2= mat.K2*Du;
       % Assume initial state of spring is elastic.
       mat.Fs1 = 0.;
       mat.Fs2 = 0.;     
       mat.YieldCode = 0;   
       mat.xprev = 0;% material initial deformation
       mat.Fs = mat.Fs1+mat.Fs2+ mat.K3*mat.xprev;
       
    case 6  % Stiffness and strength degradation model
        mat.Kip    = material(3);
        mat.Fyp     = material(4);
        mat.Fup     = material(5);
        mat.Uup     = material(6); 
        mat.Frp     = material(7); 
        mat.Urp     = material(8); 
        mat.Kin     = material(9);
        mat.Fyn     = material(10);
        mat.Fun     = material(11);
        mat.Uun     = material(12); 
        mat.Frn     = material(13); 
        mat.Urn     = material(14); 
    
        % initialize material properties
        mat.Us  = 0.0;
        mat.Uyp  = mat.Fyp/mat.Kip;
        mat.Uyn  = mat.Fyn/mat.Kin;
        mat.K1p  = mat.Kip;
        mat.K1n  = mat.Kin;
        mat.K2p  =(mat.Fup - mat.Fyp)/(mat.Uup-mat.Uyp);
        mat.K2n  =(mat.Fun - mat.Fyn)/(mat.Uun-mat.Uyn);
        mat.K3p  =(mat.Fup - mat.Frp)/(mat.Uup-mat.Urp);
        mat.K3n  =(mat.Fun - mat.Frn)/(mat.Uun-mat.Urn);
        mat.K4p  = 0;  % assume the stiffness is zero at residual branch
        mat.K4n  = 0;  % assume the stiffness is zero at residual branch
        mat.Ktan        = mat.K1p;	
        mat.StateCode   = 1;
        mat.Umax        = mat.Uyp;
        mat.Umin        = mat.Uyn;
        mat.Fmax        = mat.Fyp;
        mat.Fmin        = mat.Fyn;
        mat.Fs          = mat.K1p*mat.Us; 
        
        % Check input parameters
        if mat.Uup < mat.Uyp || mat.Uun > mat.Uyn
            msg ='Check SDegrading material input: Uu is less than Uy';
            errordlg(msg,'Input Error');
        end
        if mat.Uup > mat.Urp || mat.Uun < mat.Urn
            msg ='Check SDegrading material input: Ur is less than Uu';
            errordlg(msg,'Input Error');
        end
        if mat.Frp > mat.Fup || mat.Frn < mat.Fun
            msg ='Check SDegrading material input: Fr is greater than Fu';
            errordlg(msg,'Input Error');
        end 
        if mat.Fup < mat.Fyp || mat.Fun > mat.Fyn
            msg ='Check SDegrading material input: Fu is less than Fy';
            errordlg(msg,'Input Error');
        end 
        
    
        
    case 7  % Concrete01 material properties
        mat.fpc    = material(3);
        mat.epsc0     = material(4);
        mat.fpcu     = material(5);
        mat.epscu     = material(6);        
        mat.CminStrain = 0.0;
        mat.CendStrain = 0.0;   
        mat.Cstress = 0.0;        
        mat.Cstrain = 0.0;   
        
        if (material(3) > 0.0)
          mat.fpc = -mat.fpc;
        end     
        if (material(4) > 0.0)
          mat.epsc0 = -mat.epsc0;
        end   
        if (material(5) > 0.0)
          mat.fpcu = -mat.fpcu;
        end
        if (material(6) > 0.0)
          mat.epscu = -mat.epscu;
        end  
        
        % Initial tangent
        Ec0 = 2*mat.fpc/mat.epsc0;
        mat.Ttangent = Ec0;
        mat.CminStrain = 0.0;
        mat.CunloadSlope = Ec0;
        mat.CendStrain = 0.0;

        %State variables
        mat.Cstrain = 0.0;
        mat.Cstress = 0.0;
        mat.Ctangent = Ec0;
        
        
    case 8  % Concrete02 material properties            
        mat.fc    = material(3);
        mat.epsc0     = material(4);
        mat.fcu     = material(5);
        mat.epscu     = material(6); 
        mat.rat   = material(7);
        mat.ft   = material(8);
        mat.Ets   = material(9);
        
        mat.ecminP = 0.0;
        mat.deptP = 0.0;

        mat.eP = 2.0*mat.fc/mat.epsc0;
        mat.epsP = 0.0;
        mat.sigP = 0.0;
        mat.eps = 0.0;
        mat.sig = 0.0;
        mat.e = 2.0*mat.fc/mat.epsc0;
        
     case 9  % Steel02 material properties            
        mat.Fy    = material(3);
        mat.E0     = material(4);
        mat.b     = material(5);
        mat.R0     = material(6); 
        mat.cR1   = material(7);
        mat.cR2   = material(8);
        mat.a1   = material(9);   
        mat.a2   = material(10); 
        mat.a3   = material(11); 
        mat.a4   = material(12); 
        mat.sigini  = material(13); 
        
        mat.Esh = mat.b*mat.E0;
        mat.epsy = mat.Fy / mat.E0;
        
        mat.konP = 0;
        mat.kon = 0;
        mat.eP = mat.E0;
        mat.epsP = 0.0;
        mat.sigP = 0.0;
        mat.sig = 0.0;
        mat.eps = 0.0;
        mat.e = mat.E0;
        
        mat.epsmax = mat.Fy/mat.E0;
        mat.epsmin = -mat.epsmax;
        mat.epspl = 0.0;
        mat.epss0 = 0.0;
        mat.sigs0 = 0.0;
        mat.epsr = 0.0;
        mat.sigr = 0.0;
        
        mat.epsmaxP = mat.Fy/mat.E0;
        mat.epsminP = -mat.epsmaxP;
        mat.epsplP = 0.0;
        mat.epss0P = 0.0;
        mat.sigs0P = 0.0;
        mat.epssrP = 0.0;
        mat.sigsrP = 0.0;

        if (mat.sigini ~= 0.0) 
          mat.epsP = mat.sigini/E0;
          mat.sigP = mat.sigini;
        end
                
otherwise
        errordlg(['Unknown Material Type', num2str(mat.Type)],'Input Error');        
            
            
end