%% Return the Material State
function [ mat Et sc] = MaterialState( dep, mat, vc)
%MATERIALSTATE Summary of this function goes here
%   dep = incremental deformation at the current time step, except for
%   Bouc-Wen model
%   mat = material data structure
%   vc  = 0 for other materials
%       = current velocity for Bouc-Wen material
%
switch  mat.Type
    case 2
        mat = bilinearMat(dep, mat);
        Et = mat.Et;
        sc = mat.sc;
    case 3
        mat = hystereticMat(dep, mat);
        Et = mat.Ttangent;
        sc = mat.Cstress;
    case 4
        mat = BoucWen(dep, vc,  mat);
        Et = mat.k;
        sc = mat.Qj;
    case 5
        mat = trilinear(dep, mat);
        Et = mat.Kt;
        sc = mat.Fs;
    case 6
        mat = SDegrading(dep, mat);
        Et = mat.Ktan;
        sc = mat.Fs;      
    case 7
         mat = Concrete01(dep, mat);
         Et = mat.Ttangent;
         sc = mat.Cstress;
    case 8
         mat = Concrete02(dep, mat);
         Et = mat.e;
         sc = mat.sig;
    otherwise
        msg =['Assigned material is defined'];
        errordlg(msg,'Input Error');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mat] = bilinearMat(deps, mat )
%	[mat]=BILINEARMAT(deps,mat) For strain increment deps
%	and mat compute new properties for bilinear material.
%	mat.SY, mat.E, mat.EP are the yield stress,
%	first slop and second slope respectively.  materialProp.SC,
%	mat.SB, materialSE are the current stress, back stress, and
%	current strain, respectively. mat.ET is the current tangent modulus.

% Material properties
sigy = mat.sy;		% yield stress
E    = mat.e;		% E
Ep   = mat.ep;		% alpha*E

% State from last converged load step
sig  = mat.sc;		% stress
sigb = mat.sb;		% backstress
code = mat.cd;		% yield code
					%	=0, elastic
					%	=1, plastic

if code==0 || ( code==1 && (sig-sigb)*deps < 0 )
	
	% Strain to plastic loading from elastic state
	deps1 = ( sign(deps) * sigy + sigb - sig ) / E;

	% Elastic
	sig  = sig + E * deps;
	Et   = E;
	code = 0;
	
	% Check if plasic loading from elastic state
	if abs(deps) > abs(deps1)
		sig  = sig  + ( Ep - E ) * ( deps - deps1 );
		sigb = sigb +   Ep       * ( deps - deps1 );
		Et   = Ep;
		code = 1;
	end

else
	
	% Continue plastic loading
	sig  = sig  + Ep * deps;
	sigb = sigb + Ep * deps;
	Et   = Ep;
	code = 1;
	
end

% Return updated state. 
mat.sc = sig;
mat.sb = sigb;
mat.cd = code;
mat.es = mat.es + deps;	
mat.Et = Et;

%%%%%%%%%%%%%%%%%%%% End of bilinearMat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% source code for hysteretic material is taken from OpenSEES uniaxial
% materials
% 
% following variables are materialProp struct variables
%
% s1p, e1p = stress and strain at first point of the envelope
%           in the positive direction
% s2p, e2p = stress and strain at second point of the envelope
%           in the positive direction
% s3p, e3p = stress and strain at third point of the envelope
%           in the positive direction (optional)
% s1n, e1n = stress and strain at first point of the envelope
%           in the negative direction*
% s2n, e2n = stress and strain at second point of the envelope
%           in the negative direction*
% s3n, e3n = stress and strain at third point of the envelope 
%           in the negative direction (optional)*
% pinchX = pinching factor for strain (or deformation) during reloading
% pinchY = pinching factor for stress (or force) during reloading
% damage1= damage due to ductility: D1(mu-1)
% damage2= damage due to energy: D2(Eii/Eult)
% beta   = power used to determine the degraded unloading stiffness
%           based on ductility, mu-beta (optional, default=0.0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mat] = hystereticMat(deps, mat)
   
mat.Tstrain = deps + mat.Cstrain;

if (mat.TloadIndicator == 0 && mat.Tstrain == 0.0)
    return ;
end

mat.TrotMax = mat.CrotMax;
mat.TrotMin = mat.CrotMin;
mat.TenergyD = mat.CenergyD;
mat.TrotPu = mat.CrotPu;
mat.TrotNu = mat.CrotNu;
mat.TloadIndicator = mat.CloadIndicator;

if (mat.TloadIndicator == 0)
    if (deps < 0.0)
        mat.TloadIndicator = 2;
    else
        mat.TloadIndicator = 1;
    end
end
  
if (mat.Tstrain >= mat.CrotMax) 
    mat.TrotMax = mat.Tstrain;
    mat.Ttangent= posEnvlpTangent(mat);
    mat.Tstress = posEnvlpStress(mat,mat.Tstrain);
    mat.TloadIndicator=1;
elseif (mat.Tstrain <= mat.CrotMin) 
    mat.TrotMin = mat.Tstrain;
    mat.Ttangent= negEnvlpTangent(mat);
    mat.Tstress = negEnvlpStress(mat,mat.Tstrain);
    mat.TloadIndicator=2;
else 
    if (deps < 0.0)
        mat=negativeIncrement(mat, deps);
    elseif (deps > 0.0)
        mat=positiveIncrement(mat, deps);
    end
end
  
mat.TenergyD = mat.CenergyD + 0.5*(mat.Tstress + mat.Cstress)*deps;
   
% update the state for next time step
mat.CrotMax = mat.TrotMax;
mat.CrotMin = mat.TrotMin;
mat.CrotPu = mat.TrotPu;
mat.CrotNu = mat.TrotNu;
mat.CenergyD = mat.TenergyD;
mat.CloadIndicator = mat.TloadIndicator;

mat.Cstress = mat.Tstress;
mat.Cstrain = mat.Tstrain;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mat] = negativeIncrement(mat, deps)

kn = (mat.CrotMin/mat.rot1n)^mat.beta;
if kn < 1.0 kn =  1.0; 
else kn = 1.0/kn; end;

kp = (mat.CrotMax/mat.rot1p)^mat.beta;
if kp < 1.0 kp =  1.0;
else kp = 1.0/kp; end;

if (mat.TloadIndicator == 1) 
    mat.TloadIndicator = 2;
    if (mat.Cstress >= 0.0) 
    
        mat.TrotPu = mat.Cstrain - mat.Cstress/(mat.Eup* kp);
		energy = mat.CenergyD -0.5*mat.Cstress/(mat.Eup* kp)*mat.Cstress;
        damfc = 0.0;
        
        if (mat.CrotMax > mat.rot1p)
            damfc = mat.damfc1*(mat.CrotMax-mat.rot1p)/mat.rot1p + ...
                mat.damfc2* energy/ mat.energyA;
        end
        mat.TrotMin =  mat.CrotMin*(1.0+ damfc);		
    end
end

mat.TloadIndicator = 2;

if mat.TrotMin >= mat.rot1n
    mat.TrotMin =  mat.rot1n;
end

minmom = negEnvlpStress(mat, mat.TrotMin);
rotlim = posEnvlpRotlim(mat);
if (rotlim < mat.TrotPu)
    rotrel  = rotlim;
else
    rotrel = mat.TrotPu;
end

rotmp2 = mat.TrotMin - (1.0-mat.pinchY)*minmom/(mat.Eun*kn);
rotch  = rotrel + (rotmp2-rotrel)*mat.pinchX;     

if (mat.Tstrain > mat.TrotPu) 
	mat.Ttangent = mat.Eup*kp;
	mat.Tstress = mat.Cstress + mat.Ttangent*deps;
    if (mat.Tstress <= 0.0) 
        mat.Tstress = 0.0;
    	mat.Ttangent = mat.Eup*1.0e-9;
    end
elseif (mat.Tstrain <= mat.TrotPu && mat.Tstrain > rotch)
    
    if (mat.Tstrain >= rotrel) 
		mat.Tstress = 0.0;
		mat.Ttangent = mat.Eun*1.0e-9;
    else 
		mat.Ttangent = minmom *mat.pinchY/(rotch-rotrel);
		tmpmo1 = mat.Cstress + mat.Eun* kn*deps;
		tmpmo2 = (mat.Tstrain-rotrel)*mat.Ttangent;
        if (tmpmo1 > tmpmo2) 
			mat.Tstress = tmpmo1;
			mat.Ttangent = mat.Eun*kn;
		else
			mat.Tstress = tmpmo2;
        end
    end
else 
	mat.Ttangent = (1.0-mat.pinchY)*minmom/(mat.TrotMin-rotch);
	tmpmo1 = mat.Cstress + mat.Eun*kn*deps;
	tmpmo2 = mat.pinchY*minmom + (mat.Tstrain-rotch)*mat.Ttangent;
    if (tmpmo1 > tmpmo2)
		mat.Tstress = tmpmo1;
		mat.Ttangent = mat.Eun*kn;
    else
        mat.Tstress = tmpmo2;
    end
end
%%%%%%%%%%%%%%% end of negativeIncrement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mat] = positiveIncrement(mat, deps)

kn = (mat.CrotMin/mat.rot1n)^mat.beta;
if (kn < 1.0) kn = 1.0;
else kn = 1.0/kn; end

kp = (mat.CrotMax/mat.rot1p)^mat.beta;
if (kp < 1.0) kp = 1.0 ;
else kp =1.0/kp; end

if mat.TloadIndicator == 2     
    mat.TloadIndicator = 1;
    
    if (mat.Cstress <= 0.0)         
        mat.TrotNu = mat.Cstrain - mat.Cstress/(mat.Eun*kn);
		energy = mat.CenergyD - 0.5*mat.Cstress/(mat.Eun*kn)*mat.Cstress;
		damfc = 0.0;
        if (mat.CrotMin < mat.rot1n)
			damfc = mat.damfc2*energy/mat.energyA+ ...
                mat.damfc1*(mat.CrotMin-mat.rot1n)/mat.rot1n;
        end
        % update rotMax
		mat.TrotMax = mat.CrotMax*(1.0+damfc);
    end
end

mat.TloadIndicator = 1;

if (mat.TrotMax <= mat.rot1p) 
    mat.TrotMax = mat.rot1p;
end

maxmom = posEnvlpStress(mat, mat.TrotMax);
rotlim = negEnvlpRotlim(mat);
if (rotlim > mat.TrotNu) rotrel= rotlim;
else rotrel= mat.TrotNu; end;

rotmp2 = mat.TrotMax - (1.0-mat.pinchY)*maxmom/(mat.Eup*kp);
rotch = rotrel + (rotmp2-rotrel)*mat.pinchX; 

if (mat.Tstrain < mat.TrotNu)
	mat.Ttangent = mat.Eun*kn;
	mat.Tstress = mat.Cstress + mat.Ttangent*deps;
    if (mat.Tstress >= 0.0) 
        mat.Tstress = 0.0;
		mat.Ttangent = mat.Eun*1.0e-9;
    end
elseif (mat.Tstrain >= mat.TrotNu && mat.Tstrain < rotch) 
    if (mat.Tstrain <= rotrel) 
		mat.Tstress = 0.0;
		mat.Ttangent = mat.Eup*1.0e-9;
	else 
		mat.Ttangent = maxmom*mat.pinchY/(rotch-rotrel);
		tmpmo1 = mat.Cstress + mat.Eup*kp*deps;
		tmpmo2 = (mat.Tstrain-rotrel)*mat.Ttangent;
        if (tmpmo1 < tmpmo2) 
			mat.Tstress = tmpmo1;
			mat.Ttangent = mat.Eup*kp;
		else
			mat.Tstress = tmpmo2;
        end
    end
else
	mat.Ttangent = (1.0-mat.pinchY)*maxmom/(mat.TrotMax-rotch);
	tmpmo1 = mat.Cstress + mat.Eup*kp*deps;
	tmpmo2 = mat.pinchY*maxmom + (mat.Tstrain-rotch)*mat.Ttangent;
    if (tmpmo1 < tmpmo2) 
		mat.Tstress = tmpmo1;
		mat.Ttangent = mat.Eup*kp;
	else
		mat.Tstress = tmpmo2;
    end
end
%%%%%%%%%%%%%%%%%%%% End of positiveIncrement %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
function [posTangent] = posEnvlpTangent(mat)   
    
if (mat.Tstrain <= 0.0)
    posTangent= mat.E1p*1.0e-9;
elseif (mat.Tstrain <= mat.rot1p)
    posTangent= mat.E1p;
elseif (mat.Tstrain <= mat.rot2p)
    posTangent= mat.E2p;
elseif (mat.Tstrain <= mat.rot3p || mat.E3p > 0.0)
    posTangent= mat.E3p;
else
    posTangent= mat.E1p*1.0e-9;
end
%%%%%%%%%%%%%%%%%%%% End of posEnvlpTangent %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [negTangent] = negEnvlpTangent( mat)   

if (mat.Tstrain > 0.0)
    negTangent= mat.E1n*1.0e-9;
elseif (mat.Tstrain >= mat.rot1n)
    negTangent= mat.E1n;
elseif (mat.Tstrain >= mat.rot2n)
    negTangent= mat.E2n;
elseif (mat.Tstrain >= mat.rot3n || mat.E3n > 0.0)
    negTangent= mat.E3n;
else
    negTangent= mat.E1n*1.0e-9;
end
%%%%%%%%%%%%%%%%%%% End of negEnvlpTangent %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [posStress] = posEnvlpStress(mat, strain)     
    
if (strain <= 0.0)
    posStress = 0.0;
elseif (strain <= mat.rot1p)
    posStress = mat.E1p*strain;
elseif (strain <= mat.rot2p)
    posStress = mat.mom1p + mat.E2p*(strain-mat.rot1p);
elseif (strain <= mat.rot3p || mat.E3p > 0.0)
    posStress = mat.mom2p + mat.E3p*(strain-mat.rot2p);
else
    posStress = mat.mom3p;
end
%%%%%%%%%%%%%%%%%%%% End of posEnvlpStress %%%%%%%%%%%%%%%%%%%%%%%%
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [negStress] = negEnvlpStress(mat, strain)    
    
if (strain >= 0.0)
    negStress = 0.0;
elseif (strain >= mat.rot1n)
    negStress = mat.E1n*strain;
elseif (strain >= mat.rot2n)
    negStress = mat.mom1n + mat.E2n*(strain-mat.rot1n);
elseif (strain >= mat.rot3n || mat.E3n > 0.0)
    negStress = mat.mom2n + mat.E3n*(strain-mat.rot2n);
else
    negStress = mat.mom3n;
end
%%%%%%%%%%%%%%%%%%%% End of negEnvlpStress %%%%%%%%%%%%%%%%%%%%%%%%  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neglimit] = negEnvlpRotlim(mat)

neglimit =  -1.0e16;

if (mat.CrotMin >= mat.rot1n) 
    neglimit = -1.0e16; 
end;
if (mat.CrotMin < mat.rot1n && mat.CrotMin >= mat.rot2n && mat.E2n < 0.0)
    neglimit = mat.rot1n - mat.mom1n/mat.E2n;
end;
if (mat.CrotMin < mat.rot2n && mat.E3n < 0.0) 
    neglimit = mat.rot2n - mat.mom2n/mat.E3n;
end;

if (neglimit ==  -1.0e16) 
    neglimit = -1.0e16; 
elseif (negEnvlpStress(mat, neglimit) < 0)
    neglimit = -1.0e16; 
end
%%%%%%%%%%%%%%%%%%%% End of negEnvlpRotlim %%%%%%%%%%%%%%%%%%%%%%%%  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [poslimit] = posEnvlpRotlim(mat)

poslimit =  1.0e16;

if (mat.CrotMax <= mat.rot1p) 
    poslimit =  1.0e16; 
end;
if (mat.CrotMax > mat.rot1p && mat.CrotMax <= mat.rot2p && mat.E2p < 0.0)
    poslimit = mat.rot1p - mat.mom1p/mat.E2p;
end;
if (mat.CrotMax > mat.rot2p && mat.E3p < 0.0)
    poslimit = mat.rot2p - mat.mom2p/mat.E3p;
end;

if (poslimit == 1.0e16)
    poslimit = 1.0e16;
elseif (posEnvlpStress(mat,poslimit) > 0)
    poslimit = 1.0e16;
end
%%%%%%%%%%%%%%%%%%%% End of posEnvlpRotlim %%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%% End of hystereticMat %%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [BW]= BoucWen( dx, vcur, BW)
% ---------------------------------------------
% x,v: displacement and velocity
%   Q: restoring force
%   z: hysteretic component
%  zj: z in next step, i.e. z_i+1
%  BW: numerous parameters for Bouc-Wen model
%   i: step number

    xcur =  dx + BW.Uprev;
% Fourth order R-K forward integration
    dz1 = dx*f(BW,vcur,BW.zj);
    dz2 = dx*f(BW,vcur,BW.zj+dz1/2);
    dz3 = dx*f(BW,vcur,BW.zj+dz2/2);
    dz4 = dx*f(BW,vcur,BW.zj+dz3);
    zj = BW.zj + (dz1+2*dz2+2*dz3+dz4)/6;
 
umax=(abs(BW.Dmin)+ BW.Dmax)*0.5;    
BW.k=BW.k1*exp(-umax/BW.uref) + BW.k2;    
BW.Qj = BW.c*(abs(vcur)^BW.a)*sign(vcur) + BW.alpha*BW.k*xcur +...
    (1-BW.alpha)*BW.k*BW.uy*zj;

% update deformation for next time step
BW.Uprev = xcur;
BW.Dmax= max(xcur,BW.Dmax); 
BW.Dmin= min(xcur,BW.Dmin); 
BW.zj =  zj;

function y=f(BW,v,z)
% This function describes the nonlinear differential equation that governs
% the hysteretic component, z.
% The returned value is "dz/dx" or "z_dot/x_dot"
y =(1/BW.uy)*(1-BW.beta*sign(v)*abs(z)^BW.n*sign(z)-BW.gamma*abs(z)^BW.n);
%%%%%%%%%%%%%%%%%%% End of Bouc-Wen Model %%%%%%%%%%%%%%%%%%%%%%%%  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [pro] = trilinear(dx, pro)

% calculate current deformation 
    x = pro.xprev + dx; 

if pro.YieldCode == 0

	pro.Fs1 = pro.Fs1+pro.K1*dx;
    pro.Fs2 = pro.Fs2+pro.K2*dx;
	% if it is beyond elastic range, limit to yield force
	% and set m_nYieldCode to 1
	if pro.Fs1 > pro.Fy1
		
        pro.Fs1 = pro.Fy1 ;
		pro.YieldCode = 1;
        
        if(pro.Fs2> pro.Fy2)
			pro.Fs2 = pro.Fy2 ;
		    pro.YieldCode = 2;
        end
    elseif pro.Fs1 < -pro.Fy1
		
        pro.Fs1 = -pro.Fy1;
		pro.YieldCode = 1 ;
       
		if pro.Fs2 < -pro.Fy2
			pro.Fs2 = - pro.Fy2 ;
		    pro.YieldCode = 2;
        end
    end
% 1st post-yielding region    
elseif pro.YieldCode == 1
    
    % assume spring2 remain elastic
    pro.Fs2= pro.Fs2+pro.K2*dx;
    
    % keep loading
    if pro.Fs1*dx > 0.
        if pro.Fs2 > pro.Fy2
            pro.Fs2 = pro.Fy2;
            pro.YieldCode=2;            
        elseif pro.Fs2 < -pro.Fy2
            pro.Fs2 = -pro.Fy2;
            pro.YieldCode=2;
        end
        % unloading
    elseif pro.Fs1*dx < 0.
        pro.Fs1= pro.Fs1 + pro.K1*dx;
        pro.YieldCode = 0;
        
        if pro.Fs1 < -pro.Fy1 
            pro.Fs1 = -pro.Fy1;
            pro.YieldCode = 1;
            if pro.Fs2 < -pro.Fy2
                pro.Fs2 = -pro.Fy2;
                pro.YieldCode = 2;
            end
        elseif pro.Fs1 > pro.Fy1
            pro.Fs1 = pro.Fy1;
            pro.YieldCode = 1;
            if pro.Fs2 > pro.Fy2
                pro.Fs2 = pro.Fy2;
                pro.YieldCode = 2;
            end
        end
    end
elseif pro.YieldCode ==2
    % unloading
	if (pro.Fs1*dx) < 0.
		pro.Fs1 = pro.Fs1 + pro.K1*dx;
		pro.Fs2 = pro.Fs2 + pro.K2*dx;
		pro.YieldCode = 0;
		
        if pro.Fs1 < -pro.Fy1
            pro.Fs1 = -pro.Fy1;
			pro.YieldCode = 1;
            if pro.Fs2 < -pro.Fy2
				pro.Fs2 = -pro.Fy2;
				pro.YieldCode = 2;
            end            
        elseif pro.Fs1 > pro.Fy1
			pro.Fs1 = pro.Fy1;
			pro.YieldCode = 1;
            if pro.Fs2> pro.Fy2
				pro.Fs2 = pro.Fy2;
				pro.YieldCode = 2;
			end
        end
    end
end

pro.Fs = pro.Fs1 + pro.Fs2 + pro.K3*x;

% updata current state for next time step
if pro.YieldCode == 0
    pro.Kt = pro.K1 +pro.K2 +pro.K3;
elseif pro.YieldCode == 1
    pro.Kt = pro.K2 +pro.K3;    
else
    pro.Kt = pro.K3;
end
pro.xprev = x;
%%%%%%%%%%%%%%%%%%%%%%%%% End of Trilinear %%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [mat]= SDegrading(dx, mat)
%	[mat]=SDegrading(dx,mat) For incremental deformation, dx
%	and mat compute new properties for SDegrading material.
%	mat variables are defined in CreateMaterial.m
%   Line search algorithm is used to search for the current state

% local variables
% accFac	: accumulated factor to trace out nonlinear path
% factor	: factor used to trace out nonlinear path
% fac	    : factor used to trace out nonlinear path
% dF        :
% uZero     :
% 

dF = 0;	
accFac = 0.;

while (accFac < 0.99999999)
    
    factor = 1.0 - accFac;
    % current state is in elastic branch
    if (mat.StateCode == 1)
        % elastically loading
        if (dx>0.)
            fac = (mat.Fyp-mat.Fs)/(mat.K1p*dx);
            if(fac<factor)
                factor = fac;
                mat.StateCode = 2;
            end
        elseif (dx<0.)
            % in negative spring force
            fac = -mat.Fs/(mat.K1p*dx);
            % change in state
            if (fac < factor)
                factor = fac;
                mat.StateCode = -1;
            end
        end
        dF = mat.K1p*dx;
        
    elseif (mat.StateCode == -1)
        
        % elastically loading
        if (dx<0.)
            fac = (mat.Fyn-mat.Fs)/(mat.K1n*dx);
            if(fac<factor)
                factor = fac;
                mat.StateCode = -2;
            end
            % elastically unloading
        elseif (dx>0.)
            % in negative spring force
            fac = -mat.Fs/(mat.K1n*dx);
            % change in state
            if(fac < factor)
                factor = fac;
                mat.StateCode = 1;
            end
        end
        dF = mat.K1n*dx;
        % current state is in yielding branch
    elseif (mat.StateCode == 2)
        
        % elastically unloading along the line 5
        if(dx < 0.)
            factor =0.;
            mat.StateCode = 5;
            mat.Umax = mat.Us;
            mat.Fmax = mat.Fs;
            % continue to yield
        elseif (dx > 0.)
            fac = (mat.Uup - mat.Us) / dx;
            if(fac<factor)
                factor = fac;
                mat.StateCode  = 3;
            end
            dF = mat.K2p*dx;
       %     mat.Fmax = mat.Fs;
       %     mat.Umax = mat.Us;
        end
        
    elseif (mat.StateCode == -2)
        
        % elastically unloading along the line 5
        if (dx > 0.)
            factor =0.;
            mat.StateCode = -5;
            mat.Umin = mat.Us;
            mat.Fmin = mat.Fs;
            % continue to yield
        elseif (dx <0.)
            fac = (mat.Uun - mat.Us) / dx;
            if(fac<factor)
                factor = fac;
                mat.StateCode = -3;
            end
            dF = mat.K2n*dx;
    %        mat.Fmin = mat.Fs;
    %        mat.Umin = mat.Us;
        end
        % current state is in descending branch
    elseif ( mat.StateCode == 3)
        
       %	 elastically unloading
        if (dx < 0.)
            factor =0.;
            mat.StateCode = 5;
            % in positive spring force
            mat.Umax = mat.Us;
            mat.Fmax = mat.Fs;
            % continue to yield
        elseif (dx > 0.)
            % in positive spring force
            fac = (mat.Frp-mat.Fs) /(mat.K3p*dx);
            % change in state
            if (fac<factor)
                factor = fac;
                mat.StateCode  = 4;
            end
            dF = mat.K3p*dx;
        end
    elseif ( mat.StateCode == -3)
        
         % elastically unloading
        if (dx > 0.)
            factor =0.;
            mat.StateCode = -5;
            % in positive spring force
            mat.Umin = mat.Us;
            mat.Fmin = mat.Fs;
            % continue to yield
        elseif (dx < 0.)
            % in positive spring force
            fac = (mat.Frn-mat.Fs) /(mat.K3n*dx);
            % change in state
            if (fac<factor)
                factor = fac;
                mat.StateCode  = -4;
            end
            dF = mat.K3n*dx;
        end
    elseif (mat.StateCode == 4)
        
        % elastically unloading
        if(dx < 0.)
            factor =0.;
            mat.StateCode = 5;
            mat.Umax = mat.Us;
            mat.Fmax = mat.Fs;
        elseif (dx > 0.)
            dF = mat.K4p*dx;
        end
        
    elseif (mat.StateCode == -4)
        
        % elastically unloading
        if (dx > 0.)
            factor =0.;
            mat.StateCode = -5;
            mat.Umin = mat.Us;
            mat.Fmin = mat.Fs;
        elseif (dx < 0.)
            dF = mat.K4n*dx;
        end
        % unloading branch
    elseif (mat.StateCode == 5)

        % elastically reloading branch
        if ( dx > 0.)
            fac = (mat.Fmax - mat.Fs)/(mat.K1p*dx);
            if(fac < factor)
                % update factor
                factor = fac;
                if mat.Umax >= mat.Urp
                    mat.StateCode = 4;
                elseif mat.Umax >= mat.Uup
                    mat.StateCode = 3;                    
                elseif mat.Umax >= mat.Uyp
                    mat.StateCode = 2;
                end
            end
        elseif ( dx < 0.)
            fac =  -mat.Fs/(mat.K1p*dx);
            if(fac < factor)
                factor = fac;
                mat.StateCode = -6;
                uZero = mat.Us+factor*dx;
                mat.K6	= mat.Fmin/(mat.Umin-uZero);
            end
        end
        dF = mat.K1p*dx;
        
    elseif (mat.StateCode == -5)
        
        % elastically reloading branch
        if( dx < 0.)
            fac = (mat.Fmin - mat.Fs)/(mat.K1n*dx);
            if (fac < factor)
                % update factor
                factor = fac;
                if mat.Umin <= -mat.Urn
                    mat.StateCode = -4;
                elseif mat.Umin <= -mat.Uun
                    mat.StateCode = -3;                    
                elseif mat.Umin <= -mat.Uyn
                    mat.StateCode = -2;
                end
            end
        elseif ( dx > 0.)
            fac = -mat.Fs/(mat.K1n*dx);
            if (fac < factor)
                factor = fac;
                mat.StateCode = 6;
                uZero = mat.Us+factor*dx;
                mat.K6	= mat.Fmax/(mat.Umax-uZero);
            end
        end
        dF = mat.K1n*dx;
        
        % current state is in shooting branch
    elseif (mat.StateCode == 6)
        
        % elastically unloading
        if (dx< 0.)
            factor = 0.;
            mat.StateCode = 7;
            mat.Ftem = mat.Fs;
            mat.Utem = mat.Us;
            % continue to shoot
        elseif (dx> 0.)
            fac = (mat.Umax-mat.Us)/dx;
            if (fac<factor)
                factor = fac;
                if mat.Umax >= mat.Urp
                    mat.StateCode = 4;
                elseif mat.Umax >= mat.Uup
                    mat.StateCode = 3;                    
                elseif mat.Umax >= mat.Uyp
                    mat.StateCode = 2;
                end
            end
            dF = mat.K6*dx;
        end
        
    elseif (mat.StateCode == -6)
 
        % elastically unloading
        if(dx>0.)
            factor = 0.;
            mat.StateCode = -7;
            mat.Ftem = mat.Fs;
            mat.Utem = mat.Us;
        elseif (dx<0.)
            fac = (mat.Umin-mat.Us)/dx;
            if (fac< factor)
                factor= fac;
                if mat.Umin <= mat.Urn
                    mat.StateCode = -4;
                elseif mat.Umin <= mat.Uun
                    mat.StateCode = -3;                    
                elseif mat.Umin <= mat.Uyn
                    mat.StateCode = -2;
                end
            end
            dF = mat.K6*dx;
        end
        
        % unloading branch off shooting branch
    elseif (mat.StateCode == 7)

        % elastically reloads
        if(dx>0.)
            fac = (mat.Utem -mat.Us)/dx;
            if (fac < factor)
                factor  = fac;
                mat.StateCode  = 6;
            end
            % continue unloading elastically
        elseif (dx<0.)
            fac = -mat.Fs/(mat.K1p*dx);
            % change in state
            if (fac < factor)
                factor = fac;
                mat.StateCode =-6;
                % update stiffness in line 6
                uZero = mat.Us + factor*dx; 
                mat.K6 = mat.Fmin/(mat.Umin-uZero);
            end
        end
        dF = mat.K1p*dx;
        
    elseif (mat.StateCode == -7)

        % elastically reloads
        if(dx<0.)
            fac = (mat.Utem -mat.Us)/dx;
            if(fac < factor)
                factor  = fac;
                mat.StateCode  = -6;
            end
            % continue unloading elastically
        elseif (dx>0.)
            fac = -mat.Fs/(mat.K1n*dx);
            % change in state
            if (fac < factor)
                factor = fac;
                mat.StateCode =6;
                % update stiffness in line 6
                uZero =mat.Us + factor*dx;
                mat.K6 = mat.Fmax/(mat.Umax-uZero);
            end
        end
        dF = mat.K1n*dx;
        
    else
        % do nothing
    end
        
    % update state
    mat.Us = mat.Us + factor*dx;
    mat.Fs= mat.Fs + factor*dF;
    
    accFac = accFac+factor;
end

% update stiffness at the current state
if mat.StateCode == 1
    mat.Ktan = mat.K1p;
elseif mat.StateCode == -1
    mat.Ktan = mat.K1n;    
elseif mat.StateCode == 2
    mat.Ktan = mat.K2p;
elseif mat.StateCode == -2
    mat.Ktan = mat.K2n;    
elseif (mat.StateCode) == 3
    mat.Ktan = mat.K3p;
elseif (mat.StateCode) == -3
    mat.Ktan = mat.K3n;
elseif (mat.StateCode) == 4
    mat.Ktan = mat.K4p;
elseif (mat.StateCode) == -4
    mat.Ktan = mat.K4n;
elseif abs(mat.StateCode) == 6
    mat.Ktan = mat.K6;
end
%%%%%%%%%%%%%%%%%%%% End of SDegrading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [mat]= ModifiedBW(dx, v, mat)
% % Modified Bouc-Wen model hysteresis model
% % C. P. Heine (2001). Simulated response of degrading hysteretic joints with slack behaviour. PhD dissertation. 
% % Virginial Polytechnic Institute and State University. Blacksburg, Virginia 
% % URL: http://scholar.lib.vt.edu/theses/available/etd-08092001-100756
% % Employ the 4th order Runge-Kutta forward integration
% 
% % ---------------------------------------------
% % x,v: displacement and velocity
% %   Q: restoring force
% %   z: hysteretic component
% %  zj: z in next step, i.e. z_i+1
% %  MBW: numerous parameters for Modified Bouc-Wen model
% %   i: step number
% 
%     xcur = dx+ mat.UPrev;
%     
% % Fourth order R-K forward integration
%     dz1 = dx*G(mat, v, mat.z, mat.e);
%     de1 = dx*F(mat, mat.z);
%     dz2 = dx*G(mat, v, mat.z+dz1/2, mat.e +de1/2);
%     de2 = dx*F(mat, mat.z +dz1/2);
%     dz3 = dx*G(mat, v, mat.z+dz2/2, mat.e +de2/2);
%     de3 = dx*F(mat, mat.z +dz2/2);
%     dz4 = dx*G(mat, v, mat.z+dz3, mat.e +de3);
%     de4 = dx*F(mat, mat.z +dz3);
%     
%     zj = mat.z + (dz1+2*dz2+2*dz3+dz4)/6;
%     ej = mat.e + (de1+2*de2+2*de3+de4)/6;
%     s = Slack(mat, zj); % slack growth
% 
% BW.Qj = s*mat.alpha*mat.k*xcur  + (1-mat.alpha)*mat.k*zj;
% 
% % update state 
% mat.UPrev = xcur;
% mat.Umax= max(xcur ,mat.Umax); 
% mat.Umin= min(xcur ,mat.Umin); 
% mat.absUmax=max([(abs(mat.Umin)), mat.Umax]);
% mat.e = ej;
% mat.z = zj;
% 
% function z = G(mat, v, z, eps)
% % This function describes the nonlinear differential equation that governs
% % the hysteretic component, z.
% % The returned value is "dz/dx" or "z_dot/x_dot"
% eta  = ETA(mat, eps);
% nu   = NU(mat, eps);
% h    = H(mat, z);
% a    = A(mat, eps); 
% z  = h*(1/eta)*(a-nu*(mat.beta*sign(v)*abs(z)^mat.n*sign(z)-mat.gamma*abs(z)^mat.n));
% 
% function de_du = F(mat, z)
% de_du = (1-mat.alpha)*mat.k*z;
% 
% function h = H(mat, z)
% h = 1-mat.Xi*exp(-z^2./(mat.psi0+ mat.delta_psi*abs(mat.absUmax))^2.);
% 
% function eta =ETA(mat, eps)
% eta = 1+ mat.delta_eta*eps;
% 
% function nu = NU(mat, eps)
% nu = 1+ mat.delta_nu*eps;
% 
% function a = A(mat, eps)
% a = mat.a0 - mat.delta_a*eps;
% 
% function s = Slack(mat, z)
% s = 1 - exp(-z^2./(mat.psi0+ mat.delta_psi*abs(mat.absUmax))^2.);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mat] = Concrete01(dStrain , mat)
   
% Get the total strain 
    
    mat.Tstrain = dStrain + mat.Cstrain;
    % Get the previous states
    mat.TminStrain = mat.CminStrain;
    mat.TendStrain = mat.CendStrain;
    mat.TunloadSlope = mat.CunloadSlope;
    mat.Tstress = mat.Cstress;
    mat.Ttangent = mat.Ctangent;    
    %mat.Tstrain=mat.Cstrain;
    
            
    
    
%     DBL_EPSILON=2.2204460492503131e-16;
%     if (abs(dStrain) < DBL_EPSILON)
%     return;
%     end

    % check for a quick return
    if (mat.Tstrain > 0.0)
        mat.Tstress = 1e-10;
        mat.Ttangent = 0;%1e-10;
     end
    
    % Calculate the trial state given the change in strain
    % determineTrialState (dStrain);
    
    %[mat] = determineTrialState (mat,dStrain);   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mat.TunloadSlope = mat.CunloadSlope;  
  
    tempStress = mat.Cstress + mat.TunloadSlope*mat.Tstrain - mat.TunloadSlope*mat.Cstrain;
  
    % Material goes further into compression
    if (mat.Tstrain < mat.Cstrain)   
        mat.TminStrain = mat.CminStrain;
        mat.TendStrain = mat.CendStrain;  
    
        [mat] = reload(mat);
    
        if (tempStress > mat.Tstress) 
            mat.Tstress = tempStress;
            mat.Ttangent = mat.TunloadSlope;
        end
    
        % Material goes TOWARD tension
    elseif (tempStress <= 0.0)
        mat.Tstress = tempStress;
        mat.Ttangent = mat.TunloadSlope;
    % Made it into tension
    else
        mat.Tstress = 1e-10;
        mat.Ttangent = 0;%1e-10;
    end
      
 
    
    
    % Save states
    mat.CminStrain = mat.TminStrain;
    mat.CendStrain = mat.TendStrain;
    mat.CunloadSlope = mat.TunloadSlope;
    mat.Cstress = mat.Tstress;
    mat.Ctangent = mat.Ttangent;
    mat.Cstrain = mat.Tstrain;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mat] = reload(mat)

 if (mat.Tstrain <= mat.TminStrain) 
    
    mat.TminStrain = mat.Tstrain;
    
    % Determine point on envelope
    [mat] = envelope(mat);     
    [mat] = unload(mat);     
  elseif (mat.Tstrain <= mat.TendStrain) 
    mat.Ttangent = mat.TunloadSlope;
    mat.Tstress = mat.Ttangent*(mat.Tstrain-mat.TendStrain);
  else 
    mat.Tstress = 1e-10;
    mat.Ttangent = 0;%1e-10;
  end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mat] = envelope(mat)  

  if (mat.Tstrain > mat.epsc0) 
    eta = mat.Tstrain/mat.epsc0;
    mat.Tstress = mat.fpc*(2*eta-eta*eta);
    Ec0 = 2.0*mat.fpc/mat.epsc0;  % changed
    mat.Ttangent = Ec0*(1.0-eta); % changed
  elseif (mat.Tstrain > mat.epscu) 
    mat.Ttangent = (mat.fpc-mat.fpcu)/(mat.epsc0-mat.epscu);
    mat.Tstress = mat.fpc + mat.Ttangent*(mat.Tstrain-mat.epsc0);
  else 
    mat.Tstress = mat.fpcu;
    mat.Ttangent = 1e-10;
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
function [mat] = unload(mat)   
  DBL_EPSILON=2.2204460492503131e-16;
  
  tempStrain = mat.TminStrain;
  
  if (tempStrain < mat.epscu)
    tempStrain = mat.epscu;
  end
  eta = tempStrain/mat.epsc0;
  ratio = 0.707*(eta-2.0) + 0.834;
  if (eta < 2.0)
    ratio = 0.145*eta*eta + 0.13*eta;
  end
  mat.TendStrain = ratio*mat.epsc0;
  temp1 = mat.TminStrain - mat.TendStrain;
  Ec0 = 2.0*mat.fpc/mat.epsc0;  % changed
  temp2 = mat.Tstress/Ec0;   % changed
  
  if (temp1 > -DBL_EPSILON) % temp1 should always be negative
    mat.TunloadSlope = Ec0;  % changed
  elseif (temp1 <= temp2)
    mat.TendStrain = mat.TminStrain - temp1;
    mat.TunloadSlope = mat.Tstress/temp1;
  else 
    mat.TendStrain = mat.TminStrain - temp2;
    mat.TunloadSlope = Ec0;  % changed
  end
%   function [mat] = determineTrialState(mat,dStrain)   
%   mat.TminStrain = mat.CminStrain;
%   mat.TendStrain = mat.CendStrain;
%   mat.TunloadSlope = mat.CunloadSlope;
%   
%   tempStress = mat.Cstress + mat.TunloadSlope*dStrain;
%   
%   % Material goes further into compression
%   if (mat.Tstrain <= mat.Cstrain)
%     
%     [mat] = reload(mat);
%     
%     if (tempStress > mat.Tstress) 
%       mat.Tstress = tempStress;
%       mat.Ttangent = mat.TunloadSlope;
%     end
%   
%   
%   % Material goes TOWARD tension
%   elseif (tempStress <= 0.0) 
%     mat.Tstress = tempStress;
%     mat.Ttangent = mat.TunloadSlope;
%   
%   % Made it into tension
%   else
%     mat.Tstress = 0.0;
%     mat.Ttangent = 0.0;
%   end    


%%%%%%%%%%%%%%%%%%%%%%%%End of Concrete01 Material%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mat] = Concrete02(deps , mat)

       ec0 = mat.fc * 2 / mat.epsc0;
       
       mat.ecmin=mat.ecminP;
       mat.dept=mat.deptP;
       
       mat.eps=mat.epsP;
       mat.sig = mat.sigP;
       mat.e = mat.eP;
        
       
       mat.eps = deps + mat.epsP;
       
       DBL_EPSILON=2.2204460492503131e-16;
       if (abs(deps) < DBL_EPSILON)
        return;
       end
       
      if (mat.eps < mat.ecmin)
        [mat.sig,mat.e] = Compr_Envlp(mat.eps,mat); 
        mat.ecmin = mat.eps;
      else
        epsr = (mat.fcu - mat.rat * ec0 * mat.epscu) / (ec0 * (1.0 - mat.rat));
        sigmr = ec0 * epsr;  
                    
       [sigmm,~] = Compr_Envlp(mat.ecmin,mat); 
       
       er = (sigmm - sigmr) / (mat.ecmin - epsr);
       ept = mat.ecmin - sigmm / er;
       
     if (mat.eps <= ept)
         sigmin = sigmm + er * (mat.eps - mat.ecmin);
         sigmax = er * 0.5 * (mat.eps - ept);   
         mat.sig = mat.sigP + ec0 * deps;
         mat.e = ec0;
         if (mat.sig <= sigmin) 
              mat.sig = sigmin;
              mat.e = er;
         end
         if (mat.sig >= sigmax) 
             mat.sig = sigmax;
             mat.e = 0.5 * er;
         end
       else 
          epn = ept + mat.dept;    
        if (mat.eps <= epn) 
          [sicn,mat.e] = Tens_Envlp(mat.dept,mat);
           if (mat.dept ~= 0.0) 
                mat.e = sicn / mat.dept;
           else 
               mat.e = ec0;
           end
        mat.sig = mat.e * (mat.eps - ept);
       else
         epstmp = mat.eps - ept;
        [mat.sig,mat.e] = Tens_Envlp(epstmp,mat);
        mat.dept = mat.eps - ept;   
        end
     end
      end
     mat.ecminP=mat.ecmin;
     mat.sigP = mat.sig;
     mat.epsP = mat.eps;
     mat.deptP = mat.dept;
     mat.eP = mat.e;
    
    

  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[sigc,Ect] = Tens_Envlp (epsc,mat)
         
  Ec0  = 2.0*mat.fc/mat.epsc0;
  
  eps0 = mat.ft/Ec0;
  epsu = mat.ft*(1.0/mat.Ets+1.0/Ec0);
  if (epsc<=eps0) 
    sigc = epsc*Ec0;
    Ect  = Ec0;
  else
    if (epsc<=epsu) 
      Ect  = -mat.Ets;
      sigc = mat.ft-mat.Ets*(epsc-eps0);
    else 
      Ect  = 1.0e-10;
      sigc = 0.0;
    end
  end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[sigc,Ect] = Compr_Envlp(epsc,mat)
    Ec0  = 2.0*mat.fc/mat.epsc0;
    ratLocal = epsc/mat.epsc0;
    
    if (epsc>=mat.epsc0)
      sigc = mat.fc*ratLocal*(2.0-ratLocal);
      Ect  = Ec0*(1.0-ratLocal);
    else
      if (epsc>mat.epscu)
         sigc = (mat.fcu-mat.fc)*(epsc-mat.epsc0)/(mat.epscu-mat.epsc0)+mat.fc;
         Ect  = (mat.fcu-mat.fc)/(mat.epscu-mat.epsc0);  
      else    
         sigc = mat.fcu;
         Ect  = 1.0e-10;     
      end
   end


        
      
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
