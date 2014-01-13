%% Create Section object
function section = CreateSection(s, matid, mat)
% for WFsection geometry
%
% s(1)  = Section ID
% s(2)  = Section Type
%       = 1, WF section 
%       = 2, general section
% s(3:n)= Section geometry information
% s(n+1)= material type 
% if s(n+1) == 1 , bilinear material 
%    s(n+2:end) = E, sigmaY, second slop ratio
% if s(n+1) == 2 , hysteretic material  
%    etc., 

section.ID   = s(1);
section.Type = s(2);

switch (section.Type)
    case 1
        % Section geometry
        h  = s(3);      % height
        bf = s(4);      % flange width
        tf = s(5);      % flange thickness
        tw = s(6);      % web thickness
        % Section discretization
        nf = s(7);		% fibers in each flange
        nw = s(8);		% fibers in web
        if (length(s)>9)
            section.GA = s(10);
        end
        
        % Properties of fibers
        fiber.yi   = [ 0 0 ];
        fiber.A    =  0;
        fiber.prop = mat{1,matid};
        
        % Create flange fibers
        fibthick  = tf / nf;
        fiber.A  = fibthick * bf;
        yf = h/2 - fibthick/2;
        for i=1:2:nf*2
            % top flange
            fiber.yi = [ -yf 1 ];
            fibers(i) = fiber;
            
            % bottom flange
            fiber.yi = [  yf 1 ];
            fibers(i+1) = fiber;
            
            yf = yf - fibthick;
        end
        
        % Create web fibers
        fibthick  = ( h - 2.*tf ) / nw;
        fiber.A  = fibthick * tw;
        yf = h/2. - tf - fibthick/2.;
        for i=1:nw
            fiber.yi = [ -yf 1 ];
            fibers(i+2*nf) = fiber;
            yf = yf - fibthick;
        end
    case 2
%            Rectangular Section geometry
%                        y
%                        ^    coverz
%                        |    |    |
%             ----------------------    --    --
%             |   o     o     o    |     |    -- covery
%             |                    |     |
%             |   o           o    |     |
%      z <--- |         +          |     H
%             |   o           o    |     |
%             |                    |     |
%             |   o     o     o    |     |    -- cover
%             ----------------------    --    --
%             |-------- B --------|
%         
        H  = s(3);      % height
        B  = s(4);      % width
        covery = s(5);   % cover thickness in y-direction
        coverz = s(6);   % cover thickness in z-direction
        coverY = H/2.0;  % The distance from the section z-axis to the edge of the cover concrete 
        coverZ = B/2.0;  % The distance from the section y-axis to the edge of the cover concrete
        coreY  = coverY - covery;
        coreZ  = coverZ - coverz;
        %confZ  = coverZ - coverz; % The distance from the section y-axis to the edge of the core concrete
        
        % Section discretization
        nfcorey = s(7);		% fibers in y direction
        nfcorez = s(8);		% fibers in z direction       
        nfcovery = s(9);	% fibers in y direction
        nfcoverz = s(10);	% fibers in z direction    
        
        % Reinforced bar 
        nLayers = s(11);		% fibers in y direction
            
        m=1;
        for k=1:nLayers % assume more than one RB layer
            area    = s(10+ 2*k);  % section area
            nBars   = s(11+ 2*k);  % number of bars
            
            for l =1:nBars
                Pos(m,1) = coreY*(1 - 2*(k-1)/(nLayers-1));              
                Pos(m,2) = coreZ*(1 - 2*(l-1)/(nBars-1));
                As(m,1)  = area;
                m = m + 1;
            end            
        end
                
        fibers = CellPatch(As, Pos, mat{1,matid(3)});
        
        coord = [-coreY -coreZ; coreY -coreZ; coreY coreZ; -coreY coreZ]; 
        fibers = horzcat(fibers, QuadPatch(nfcorey, nfcorez, coord, mat{1,matid(1)}));  
        
        coord = [coreY -coreZ; coverY -coverZ; coverY coverZ; coreY coreZ];         
        fibers = horzcat(fibers, QuadPatch(nfcovery, nfcorez, coord, mat{1,matid(2)})); 
        
        coord = [-coverY -coverZ; coverY -coverZ; coreY -coreZ; -coreY -coreZ];         
        fibers = horzcat(fibers, QuadPatch(nfcorey, nfcoverz, coord, mat{1,matid(2)})); 
        
        coord = [-coverY -coverZ; -coreY -coreZ; -coreY coreZ; -coverY coverZ];         
        fibers = horzcat(fibers, QuadPatch(nfcovery, nfcorez, coord, mat{1,matid(2)})); 
        
        coord = [-coreY coreZ; coreY coreZ; coverY coverZ; -coverY coverZ];         
        fibers = horzcat(fibers, QuadPatch(nfcorey, nfcoverz, coord, mat{1,matid(2)}));  
    
        %layer straight IDreinf barAreaSec coreY confZ coreY -confZ
end   

% Create section object
section.fibers = fibers;	% Fiber data
section.numfibers = length(fibers);
section.vs	   = [0 0]';	% Section deformations
end

function fibers = QuadPatch(nDivIJ, nDivJK, vertCoord, mat)

deltaXi  = 2.0 / nDivIJ;
deltaEta = 2.0 / nDivJK;
cVertCd  = zeros(4,2);

k=1;
for j= 1:nDivJK
    
    for i = 1:nDivIJ
        % compute natural coordinates
        
        cVertCd(1,1) = -1.0 + deltaXi  * (i-1);
        cVertCd(1,2) = -1.0 + deltaEta * (j-1);
        cVertCd(2,1) = -1.0 + deltaXi  * i;
        cVertCd(2,2) = cVertCd(1,2);
        cVertCd(3,1) = cVertCd(2,1);
        cVertCd(3,2) = -1.0 + deltaEta * j;
        cVertCd(4,1) = cVertCd(1,1);
        cVertCd(4,2) = cVertCd(3,2);
        
        % map to cartesian coordinates using bilinear shape functions
        for r = 1: 4
            
            xi  = cVertCd(r,1);
            eta = cVertCd(r,2);
            N(1) = (1.0 - xi)*(1.0 - eta)/4.0;
            N(2) = (1.0 + xi)*(1.0 - eta)/4.0;
            N(3) = (1.0 + xi)*(1.0 + eta)/4.0;
            N(4) = (1.0 - xi)*(1.0 + eta)/4.0;
            cVertCd(r,1) = 0.0;
            cVertCd(r,2) = 0.0;
            
            for s= 1:4
                cVertCd(r,1) =cVertCd(r,1) + N(s)*vertCoord(s,1);
                cVertCd(r,2) =cVertCd(r,2) + N(s)*vertCoord(s,2);
            end
            
        end
        
        area = ((cVertCd(3,1)-cVertCd(2,1))*(cVertCd(1,2)-cVertCd(2,2))...
            - (cVertCd(1,1)-cVertCd(2,1))*(cVertCd(3,2)-cVertCd(2,2))...
            + (cVertCd(1,1)-cVertCd(4,1))*(cVertCd(3,2)-cVertCd(4,2))...
            - (cVertCd(3,1)-cVertCd(4,1))*(cVertCd(1,2)-cVertCd(4,2)))/2.0; 
        cntr=[0 ;0];
        
        for m =1:4            
            m1 = rem(m,4)+1;
            ym  = cVertCd(m,1);
            zm  = cVertCd(m,2);
            ym1 = cVertCd(m1,1);
            zm1 = cVertCd(m1,2);

            dym = ym1 - ym;
            dzm = zm1 - zm;   
            integ = ym*zm + (ym*dzm + zm*dym)/2.0 + dym*dzm/3.0;

            cntr(1,1) = cntr(1,1) - dym * integ;
            cntr(2,1) = cntr(2,1) + dzm * integ;
        end
        
        cntr = cntr/area;
        fiber.yi = [cntr(1,1) 1];
        fiber.zi = [cntr(2,1) 1];
        fiber.prop = mat;
        fiber.A = area;
        
        fibers(k) = fiber;
        k=k+1;
    end
end
end

function cells = CellPatch(As, coord, mat)
    
    for m=1:length(As)
        cell.yi = [ coord(m,1) 1];
        cell.zi = [ coord(m,2) 1];
        cell.prop = mat;
        cell.A = As(m,1);        
        cells(m) = cell;
    end
end
  
      



