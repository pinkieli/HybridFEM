%% Node object is created with coordinates
function vars = CreateNode(ID, Xcoord, Ycoord, Zcoord)

% Set node coordinates
vars.ID = ID;
vars.Xcoord = Xcoord;
vars.Ycoord = Ycoord;
vars.Zcoord = Zcoord;
%vars.UX=UX; % no longer needed as of Oct. 2010, since V 4.2.3
%vars.UY=UY;
%vars.THETA=THETA;