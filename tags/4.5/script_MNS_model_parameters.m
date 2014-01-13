%% Parameters for MNS MR damper model
%
% Target MR damper: Large-scale MR damper manufactured by Lord Corp.
% (Part #: RD-9001-01, manufacturing date: Jul 15, 2005)
% Damper 1 with current input from AMC device
%
% Unit (See Table 1. in Ref.)
% c: kN*sec/m
% k: kN/m
% vm: kN*sec^2/m (corresponding to m0 in Ref.)
% vcr1, vcr2: m/sec (corresponding to xd_t+ or xd_t- in Ref.)
% a1, a4: kN  (corresponding to 'a' in Ref.)
% a2, a5: kN*sec/m (corresponding to 'b' in Ref.)
% a3, a6: non-dimensional (corresponding to 'n' in Ref.)
%
% a1, a2, a3, vcr1 ==> for positive force post-yield curve
% a4, a5, a6, vcr2 ==> for negative force post-yield curve
%
% Ref.: DEVELOPMENT OF A LARGE-SCALE MR DAMPER MODEL FOR SEISMIC
%       HAZARD MITIGATION ASSESSMENT OF STRUCTURES
%       by Yunbyeong Chae, James M. Ricles and Richard Sause
%       Proceedings of the 9th US National and 10th Canadian Conference on
%       Earthquake Engineering, Toronto, Canada, 2010.
%
% Written by Yunbyeong Chae, Mar. 2010


%% 0.0A 
c=10000; k=100000; vm=0.5; vcr1=0.01; vcr2=-0.01; 
a1=7.5; a2=243.457; a3=1.6207;    a4=-7.3; a5=-235.6136; a6=1.5979; 

PAR(1,:)=[0 c k vm vcr1 vcr2 a1 a2 a3 a4 a5 a6];  % the first column means corresponding current for parameters

%% 0.5A 
c=11000; k=100000; vm=0.5; vcr1=0.01; vcr2=-0.01; 
a1=53.1385; a2=162.5487; a3=0.8500;    a4=-a1; a5=-a2; a6=a3; 

PAR(2,:)=[0.5 c k vm vcr1 vcr2 a1 a2 a3 a4 a5 a6];

%% 1.0A 
c=12000; k=118000; vm=1.6; vcr1=0.01; vcr2=-0.01; 
a1=91.4577; a2=122.5271; a3=0.5172;    a4=-95.9947; a5=-134.8932; a6=0.6019; 

PAR(3,:)=[1 c k vm vcr1 vcr2 a1 a2 a3 a4 a5 a6];

%% 1.5A 
c=12000; k=118000; vm=1.5; vcr1=0.01; vcr2=-0.01; 
a1=126.6546; a2=152.0959; a3=0.5814;    a4=-a1; a5=-a2; a6=a3; 

PAR(4,:)=[1.5 c k vm vcr1 vcr2 a1 a2 a3 a4 a5 a6];

%% 2.0A 
c=11491; k=110030; vm=1.0489; vcr1=0.0028; vcr2=-0.0029; 
a1=148.5253; a2=166.2735; a3=0.6558;    a4=-146.7765; a5=-182.1075; a6=0.7084; 

PAR(5,:)=[2 c k vm vcr1 vcr2 a1 a2 a3 a4 a5 a6];

%% 2.5A 
c=12278; k=112890; vm=1.0418; vcr1=0.017; vcr2=-0.012;
a1=138.5239; a2=161.7807; a3=0.4564;    a4=-133.5350; a5=-171.8329; a6=0.4552; 

PAR(6,:)=[2.5 c k vm vcr1 vcr2 a1 a2 a3 a4 a5 a6];

%% 20.0A : fictitious
c=12511; k=115060; vm=1.0922; vcr1=0.01; vcr2=-0.0114; 
a1=148.0012; a2=172.0910; a3=0.5826;    a4=-133.0738; a5=-172.6988; a6=0.4437; 

PAR(7,:)=[20 c k vm vcr1 vcr2 a1 a2 a3 a4 a5 a6];

MNS.PAR = PAR;

%% Parameters for dynamics of MR damper (current: A)
MNS.alpha0=24.96;  
MNS.alpha1=3.57; 
MNS.a_pos = 0.31; 
MNS.a_neg=0.30; 

%% other parameters
MNS.TOL=1.5;  % tolerence for checking mode change (kN): default=1.5

% Initial conditions 
MNS.y0 = 0; 
MNS.yd0=0; 
MNS.z0= MNS.y0; 
MNS.zd0=0;
