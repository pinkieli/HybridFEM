%% Load Input File which has descriptions of the nodes, elements, materials
%  and integrator
function [Model, Elements, Integrator, Materials, Sections, Nodes, choice] = LoadConfiguration(INP_FILE, check)
FileToOpen = INP_FILE;              % input file path
InputFileN  = fopen(FileToOpen,'r');

%% Initialize variables
NN = 0; NM = 0; NE = 0; NSEC = 0; ND = 0; NSLAVED = 0; NGN = 0; 
lineNo = 0; iCons = 0; Model.RigidLinkNo = 0; iiNode = 0; Model.RecordElements = 'None';
format = ' '; choice = ' ';

%% Read each line of input file and purse the data till the end of the file
tline = fgetl(InputFileN); lineNo = lineNo + 1;
while ischar(tline)
    [head, remain] = strtok(tline);
    switch head
        case '#'
            % do nothing for comment lines
            
        case 'node'
            NN = NN + 1; % counts node numbers
            data = str2num(remain);
            Nodes(NN) = CreateNode(data(1), data(2), data(3), data(4));
            
        case 'fix'
            data = str2num(remain);
            Model.BOUND(data(1),:)= data(2:end);
            ND = ND + nnz(data(2:end)); % counts restrained DOFs numbers
            
        case 'equalDOF'
            iCons = iCons + 1; % counts number of equalDOF lines
            data = str2num(remain);
            NSLAVED = NSLAVED + nnz(data(3:end));
            constraintData(iCons,:) = data;
            
        case 'rigidLink'
            iCons = iCons + 1; % counts number of equalDOF lines
            Model.RigidLinkNo = Model.RigidLinkNo + 1; % counts number of rigidLink lines
            [LinkType, remain] = strtok(remain);
            switch LinkType
                case 'beam'
                    data = str2num(remain);
                    NSLAVED = NSLAVED + 3;
                    constraintData(iCons,:) = [data 1 1 1];
                    Model.RigidLinkNodeID(Model.RigidLinkNo,:)=data; % [master node #, slave node #]
                    Model.RigidLinkMaster{Model.RigidLinkNo}=Nodes(data(1));
                    Model.RigidLinkSlave{Model.RigidLinkNo}=Nodes(data(2));
                case 'bar'
                    % to be added in future
                    errordlg('Not supported in this version of HFEM')
                    choice = 'No';
                    break
                otherwise
                    errordlg(['Specify rigid link type in line no. ',num2str(lineNo)])
                    choice = 'No';
                    break
            end
            
        case 'uniaxialMaterial'
            NM = NM + 1; % count material numbers
            data = remain;
            Materials{NM} = CreateMaterial(data);
            
        case 'section'
            NSEC = NSEC + 1; % count section numbers
            data = str2num(remain);
            [par, matID] = SectionInfo(data, Materials);
            % Create a new section object
            Sections(NSEC) = CreateSection(par,matID, Materials);
            
        case 'Gsection'
            nf = 0;
            NSEC = NSEC + 1; % count section numbers
            data = str2num(remain);
            fiberDataVec = fscanf(InputFileN,'%f%f%f');
            fiberDataMat = reshape(fiberDataVec,3,length(fiberDataVec)/3)';
            lineNo = lineNo + size(fiberDataMat,1);
            Sections(NSEC) = CreateGSection(data(1), data(2), fiberDataMat, Materials, check);
            clear fiberDataVec fiberDataMat
            
        case 'element'
            NE = NE + 1; % count element numbers
            ElementData{NE} = remain;
            
        case 'gravityLoad'
            NGN = NGN + 1; % count number of gravity nodes
            gravityData(NGN,:) = str2num(remain);
            
        case 'rayleigh'
            data = str2num(remain);
            Model.T1 = data(1);
            Model.T2 = data(2);
            Model.DampingRatio = data(3);
            
        case 'interfaceNode'
            iiNode = iiNode + 1;
            intfNodeMatrix(iiNode,:) = str2num(remain);
            format = strcat(format, '%f');
            
        case 'CexpRayleigh'
            data = fscanf(InputFileN,format); lineNo = lineNo + iiNode;
            switch strtok(remain)
                case 'full' % for full matrix
                    Model.CexpRayleigh = reshape(data, iiNode, iiNode)';
                case 'lowerTriangular' % for symmetric matrix
                    cum_i = 0;
                    for i = 1 : iiNode
                        cum_i = cum_i + i;
                        Model.CexpRayleigh(i,1:i) = data(cum_i - i +1 : cum_i); 
                    end
                    for i = 1 : (iiNode - 1)
                        for j = (i + 1) : iiNode
                            Model.CexpRayleigh(i,j) = Model.CexpRayleigh(j,i);
                        end
                    end
                otherwise
                    errordlg(['Syntax error in line no. ',num2str(lineNo)])
                    choice = 'No';
                    break
            end
            
        case 'CexpInitial'
            data = fscanf(InputFileN,format); lineNo = lineNo + iiNode;
            switch strtok(remain)
                case 'full' % for full matrix
                    Model.CexpInitial = reshape(data, iiNode, iiNode)';
                case 'lowerTriangular' % for symmetric matrix
                    cum_i = 0;
                    for i = 1 : iiNode
                        cum_i = cum_i + i;
                        Model.CexpInitial(i,1:i) = data(cum_i - i +1 : cum_i); 
                    end
                    for i = 1 : (iiNode - 1)
                        for j = (i + 1) : iiNode
                            Model.CexpInitial(i,j) = Model.CexpInitial(j,i);
                        end
                    end
                otherwise
                    errordlg(['Syntax error in line no. ',num2str(lineNo)])
                    choice = 'No';
                    break
            end
            
        case 'KexpInitial'
            data = fscanf(InputFileN,format); lineNo = lineNo + iiNode;
            switch strtok(remain)
                case 'full' % for full matrix
                    Model.KexpInitial = reshape(data, iiNode, iiNode)';
                case 'lowerTriangular' % for symmetric matrix
                    cum_i = 0;
                    for i = 1 : iiNode
                        cum_i = cum_i + i;
                        Model.KexpInitial(i,1:i) = data(cum_i - i +1 : cum_i); 
                    end
                    for i = 1 : (iiNode - 1)
                        for j = (i + 1) : iiNode
                            Model.KexpInitial(i,j) = Model.KexpInitial(j,i);
                        end
                    end
                otherwise
                    errordlg(['Syntax error in line no. ',num2str(lineNo)])
                    choice = 'No';
                    break
            end
            
        case 'integrator'
            data = str2num(remain);
            Integrator = CreateIntegrator(data(3), data(4), data(5), data(1), data(2));
            
        case 'recorder'
            [RecorderType, remain] = strtok(remain);
            switch RecorderType
                case 'Element'
                    data = str2num(remain);
                    Model.RecordElements = data;
                case 'Node' % to be added in future
            end
            
        otherwise
            errordlg(['Syntax error in line no. ',num2str(lineNo)])
            choice = 'No';
            break
    end
    tline = fgetl(InputFileN);  lineNo = lineNo + 1;
end

%% Create structure
Model = CreateStructure(NN,NE,NM,NSEC,2,2,3,ND,NSLAVED,NGN,Model);
%% Purse the constraint DOF information
if iCons ~= 0 % if constraint present
    constraintData = sortrows(constraintData);    % sort array data
    i = 1;
    for j = 1 : size(constraintData, 1)
        for k = 1  : Model.NumDOF
            if ~ (constraintData(j, k + 2) == 0)
                if (Model.EqDOF(constraintData(j,1),k) == 0)
                    Model.EqDOF(constraintData(j,1),k)= i;
                    i=i+1;
                end
                Model.EqDOF(constraintData(j,2),k)= Model.EqDOF(constraintData(j,1),k);
            end
        end
    end
    clear constraintData
end
%% Convert RigidLinkMaster and RigidLinkSlave to Matricies for SIMULINK
if (Model.RigidLinkNo ~= 0)
    for j=1:Model.RigidLinkNo
        Model.RigidLinkMasterMatrix(j,1) = Model.RigidLinkMaster{j}.ID;
        Model.RigidLinkMasterMatrix(j,2) = Model.RigidLinkMaster{j}.Xcoord;
        Model.RigidLinkMasterMatrix(j,3) = Model.RigidLinkMaster{j}.Ycoord;
        Model.RigidLinkMasterMatrix(j,4) = Model.RigidLinkMaster{j}.Zcoord;
        Model.RigidLinkSlaveMatrix(j,1) = Model.RigidLinkSlave{j}.ID;
        Model.RigidLinkSlaveMatrix(j,2) = Model.RigidLinkSlave{j}.Xcoord;
        Model.RigidLinkSlaveMatrix(j,3) = Model.RigidLinkSlave{j}.Ycoord;
        Model.RigidLinkSlaveMatrix(j,4) = Model.RigidLinkSlave{j}.Zcoord;
    end
else
    % Need to define unused variables here if rigidLinks are absent
    Model.RigidLinkNodeID=[0 0]; % Define the variable with invalid data
    Model.RigidLinkMaster=Nodes(1); % Unused with invalid Node ID, Dummy data
    Model.RigidLinkSlave=Nodes(1); % Unused with invalid Node ID, Dummy data
    Model.RigidLinkMasterMatrix = zeros(3,3);
    Model.RigidLinkSlaveMatrix = zeros(3,3);
end
%% Number the dofs in the model considering restrained and constrained dofs
% number of slave dofs
numSDOF = nnz(Model.EqDOF) - max(max(Model.EqDOF));
% set dof_numberer
[Model, Nodes] = DOFNumberer(Model, numSDOF, Nodes);
%% Define if no "section" present
if NSEC ==0
    Sections = 0;
end
%% Parse the element information and create Element data structure
for i = 1:Model.NumElements
    data = str2num(ElementData{i});
    % Select the proper element type to configure
    switch data(2) 
        case 1  % Type 1: Elastic beam-column element
            Elements(i) = {CreateElement_Type1(data(1),data(2),...
                Nodes(data(3)), Nodes(data(4)), data(5), data(6),...
                Materials{data(7)}, data(8),data(9),data(10))};
        case 2  % Type 2: Experimental element with initial stiffness for simulation mode
            Elements(i) = {CreateElement_Type2(data(1),data(2),...
                Nodes(data(3)), Nodes(data(4)), data(5),data(6),...
                data(7),data(8))};
        case 3  % Type 3: Inelastic beam-column element with point-hinges DRAIN2DX type 2            
            Elements(i) = {CreateElement_Type3(data(1),data(2),...
                Nodes(data(3)), Nodes(data(4)), data(5), data(6),...
                Materials{data(7)}, data(8), data(9), data(10),...
                data(11), data(12), data(13),data(14),data(15))};
        case 4  % Type 4: Dummy load gravity column element
            Elements(i) = {CreateElement_Type4(data(1),data(2),...
                Nodes(data(3)), Nodes(data(4)), data(5), data(6),...
                data(7), data(8), data(9), data(10))};
		case 5  % Type 5: Rotational element with stiffness and strength deterioration
            Elements(i) = {CreateElement_Type5(data(1),data(2),...
                Nodes(data(3)), Nodes(data(4)), data(5), data(6),...
                data(7), data(8), data(9), data(10), data(11),...
                data(12), data(13),data(14),data(15), data(16),...
                data(17), data(18), data(19), data(20), data(21),...
                data(22), data(23),data(24),data(25), data(26),...
                data(27), data(28), data(29))};			
        case 6  % Type 6: Displacement-Based Inelastic fiber beam-column element             
            Elements(i) = {CreateElement_Type6(data(1),data(2),...
                Nodes(data(3)), Nodes(data(4)), Materials, data(5), data(6),...
                Sections(data(7)), data(8), data(9))};		
   		case 7  % Type 7: Force-Based Inelastic fiber beam-column element             
            Elements(i) = {CreateElement_Type7(data(1),data(2),...
                Nodes(data(3)), Nodes(data(4)), Materials, data(5), data(6),...
                Sections(data(7)), data(8), data(9), data(10),...
                data(11),data(12),data(13))};	
        case 8  % Type 8: zero length element
            Elements(i) = {CreateElement_Type8(data(1),data(2),...
                Nodes(data(3)), Nodes(data(4)), data(5), data(6),...
                data(7),Materials{data(8)})};
        case 9  % Type 9: panel zone element                 
            Elements(i) = {CreateElement_Type9(data(1),data(2),...
                Nodes(data(3)), Nodes(data(4)), Nodes(data(5)),...
                Nodes(data(6)), data(7), data(8), Materials{data(9)},...
                data(10), data(11), data(12), data(13), data(14),...
                data(15), data(16), Materials{data(17)}, data(18))};
        otherwise % Unknown case
            errordlg(['Unknown Element Type', num2str(data(2))],'Input Error');
            choice = 'No';
    end
end
clear ElementData
%% Create static loading
if NGN ~=0
    Model = CreateStaticLoading(Model, gravityData);
    clear gravityData
end
%% Create vector of interface DOFs
if iiNode ~= 0
    Model.InterfaceDOF = zeros(iiNode, 1);
    for i = 1 : iiNode
        Model.InterfaceDOF(i) = Model.DOF(intfNodeMatrix(i,2),intfNodeMatrix(i,1));
        if Model.InterfaceDOF(i) == -1
            errordlg('Constrained DOF cannot be an interface DOF');
            choice = 'No';
            break
        end
    end
    clear intfNodeMatrix
end
%% temporary code for checking model configuration visually
if strcmp(choice, 'No') == 0
% plot the model based on the information 
k=0; j=0;
for i=1:Model.NumElements

    if Elements{i}.Type == 9
        k=k+1;
        iN= Elements{i}.Nodes(1).ID;
        jN= Elements{i}.Nodes(2).ID;  
        kN= Elements{i}.Nodes(3).ID;
        lN= Elements{i}.Nodes(4).ID;
        PZCON{k,1} = [iN jN kN lN];
    else
        j=j+1;
        iN= Elements{i}.Nodes(1).ID;
        jN= Elements{i}.Nodes(2).ID;    
        CON{j,1} = [iN jN];
    end
end

ElmNd  = cat(1,CON{1:j,1});
% if k < double(1.)
%     PZNd   = cat(1,PZCON{1:k,1});
% end

for i=1:j
    X(1,i) = Nodes(ElmNd(i,1)).Xcoord;
    X(2,i) = Nodes(ElmNd(i,2)).Xcoord;
    Y(1,i) = Nodes(ElmNd(i,1)).Ycoord;
    Y(2,i) = Nodes(ElmNd(i,2)).Ycoord;
end

MarkerSize = 25;
figure (1)
line (X,Y,'Color',[0,0,1],'LineWidth',[0.5]);
axis equal;
hold on

for i=1:Model.NumNodes
    NdX(1,i) = Nodes(i).Xcoord;
    NdY(1,i) = Nodes(i).Ycoord;
end
scatter (NdX,NdY,'s' );

% Construct a questdlg with three options
if (check)
    choice = questdlg('Is the model configuration correct?','Model Checking', ...
          'No', 'Yes','Yes');
else
    choice = 'Yes';
end
if strcmp(choice, 'No')
    error('[FEM] Incorrect model');
end
% The following is added by C. Kolay to view node numbers
% Construct a questdlg with three options
if (check)
    DispNodeNum = questdlg('Display node numbers?','Node display', ...
        'No', 'Yes','Yes');
    if strcmp(DispNodeNum,'Yes')
%         hold on
        for i=1:Model.NumNodes
            NdID = Nodes(i).ID;
            text(NdX(1,i),NdY(1,i),num2str(NdID),'VerticalAlignment','bottom',...
                'HorizontalAlignment','left'); hold on
        end
    end
end
end

end

  
function [par, matID] = SectionInfo(data, Materials)

% search for the Material ID assigned to the section from material list
secType= data(2);
switch secType
    case {1} % wide flange section and General section
        matID=1;
        while 1
            if data(end) == Materials{matID}.ID
                par = data(1:end-1); 
                break;
            end
            matID=matID+1;
        end
    case 2  % Rectangular section
        par = data(1:end-3);
        
        for k=1:3
            matID(k)=1;
            while 1 %
                if data(end-3+k) == Materials{matID(k)}.ID
                    break;
                end
                matID(k)=matID(k)+1;
            end
        end
    otherwise
        error(['Invalid Material defined in section']);      
end

end


