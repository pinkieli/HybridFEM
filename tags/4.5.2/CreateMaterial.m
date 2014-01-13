%% Create Material object
% data(1) = material ID
% data(2) = material type
%         = 1 : elastic 
%         = 2 : bilinear elasto plastic
%         = 3 : hysteretic material
% data(3:n) =  material properties 
function mat = CreateMaterial(data)

% Material structs need to be predefined with empty variables for code
% generation
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
mat.uref = 0;
mat.Dmax = 0;   
mat.Dmin = 0; 
mat.zj = 0;
mat.Uprev =0;
mat.k = 0;
mat.Dmax = 0;
mat.Dmin = 0;     
mat.Qj = 0;

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
        mat.k2 = material(5);   
        mat.alpha = material(6);  
        mat.uy = material(7);    
        mat.a = material(8);
        mat.uref = 0.8*mat.uy;
        mat.beta = material(9);
        mat.gamma1 = material(10);
        mat.n = material(11);   
        mat.uref = material(12);
        mat.Dmax = 0.;   
        mat.Dmin = 0.; 
        mat.zj = 0.;
        mat.Uprev =0.;
        mat.k = material(13);  
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
       %Ku    = 0.04*Kinit;   
       Ku    = 0.00*Kinit;   
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
        mat.Us  = 0.;
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
        
    otherwise
        errordlg(['Unknown Material Type', num2str(mat.Type)],'Input Error');
end