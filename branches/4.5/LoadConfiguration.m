%% Load Input File which has descriptions of the nodes, elements, materials
%  and integrator
function [Model, Elements, Integrator, Materials, Sections, Nodes, choice] = LoadConfiguration(INP_FILE, check)
FileToOpen = INP_FILE;              % input file path
InputFile  = fopen(FileToOpen,'r');
[~] = fgets(InputFile);                   % skip 3 lines
[~] = fgets(InputFile);  
[~] = fgets(InputFile);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the Model information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = str2num(fgets(InputFile));   % get a configuration line
Model = CreateStructure(data(1),data(2),data(3),data(4),data(5),...
    data(6),data(7),data(8),data(9),data(10));
[~] = fgets(InputFile);                   % skip 2 lines
[~] = fgets(InputFile);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the nodal coordinate information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numNodes = 0;
data = fgets(InputFile);
% Create a new node with X,Y,Z coordinates 
while ~isstrprop(data, 'alpha')   
    data = str2num(data);
    numNodes = numNodes + 1;
    Nodes(numNodes) = CreateNode(data(1), data(2), data(3), data(4));    
    data = fgets(InputFile);
end
fprintf('numNodes: %i\n',numNodes);
[~] = fgets(InputFile); % "Boundary Condition" has been read, so ignore next line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the boundary condition and save it to Model data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
     data = fscanf(InputFile, '%d %d %d %d/n' );
     if isempty(data)
         break;
     end
     Model.BOUND(data(1),:)= data(2:Model.NumDOF+1)';
end    
[~] = fgets(InputFile);                   % skip 2 lines
[~] = fgets(InputFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the constraint dof information and save it to Model data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model.RigidLinkNo = 0;
Model.RigidLinkNodeID=[0 0]; % Define the variable with invalid data
Model.RigidLinkMaster{1}=Nodes(1); % Unused with invalid Node ID, Dummy data
Model.RigidLinkSlave{1}=Nodes(1); % Unused with invalid Node ID, Dummy data
if (Model.NumSlavedDOF ~= 0)
    i=1; ii=1;
    while 1        
        linedata=fgets(InputFile);
        [head, rem]=strtok(linedata);
    
        if( strcmp(head,'RigidLink') )
            tmp=str2num(rem);
            Model.RigidLinkNodeID(ii,:)=tmp; % [master node #, slave node #]
            Model.RigidLinkMaster{ii}=Nodes(tmp(1)); 
            Model.RigidLinkSlave{ii}=Nodes(tmp(2));  
            ii = ii + 1;
            tmp = [tmp 1 1 1];            
        else
            tmp = str2num(linedata);            
            if isempty(tmp)
                break;
            end
        end
        
        data(i,:) = tmp;
        i = i+1;
    end
    Model.RigidLinkNo = ii-1;
    data= sortrows(data);    % sort array data
    i=1;
    for j=1:size(data,1)
        for k=1:Model.NumDOF
            if ~ (data(j,k+2) == 0)
                if (Model.EqDOF(data(j,1),k) == 0)
                    Model.EqDOF(data(j,1),k)= i;
                    i=i+1;
                end
                Model.EqDOF(data(j,2),k)= Model.EqDOF(data(j,1),k);
            end
        end
    end
    [~] = fgets(InputFile);
else
    Model.RigidLinkNodeID=[0 0]; % Define the variable with invalid data
    Model.RigidLinkMaster=Nodes(1); % Unused with invalid Node ID, Dummy data
    Model.RigidLinkSlave=Nodes(1); % Unused with invalid Node ID, Dummy data
end
% Convert RigidLinkMaster and RigidLinkSlave to Matricies for SIMULINK
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
    % Need to define unused variables here
    Model.RigidLinkMasterMatrix = zeros(3,3);
    Model.RigidLinkSlaveMatrix = zeros(3,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number the dofs in the model considering restrained and constrained dofs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of slave dofs
numSDOF = nnz(Model.EqDOF) - max(max(Model.EqDOF));
% set dof_numberer
[Model, Nodes] = DOFNumberer(Model, numSDOF, Nodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the material information and create Material data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numMaterials = 0;
data = fgets(InputFile);
while ~isstrprop(data, 'alpha')       
    % Create a new material object   
    numMaterials = numMaterials + 1;
    Materials{numMaterials} = CreateMaterial(data);
    data = fgets(InputFile);
end
fprintf('numMaterials: %i\n',numMaterials);
[~] = fgets(InputFile); % skip next line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the section information and create Section data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numSections = 0;
data = fgets(InputFile);
while ~isstrprop(data, 'alpha') 
    data = str2num(data);  
    [par matID] = SectionInfo(data, Materials);
    % Create a new section object
    Sections(numSections) = CreateSection(par,matID, Materials);
    data = fgets(InputFile);
end
Sections = numSections;  % Can be 0 so elements know there are no sections 
fprintf('numSections: %i\n',numSections);
[~] = fgets(InputFile); % skip next line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the element information and create Element data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Model.NumElements
    data = str2num(fgets(InputFile));
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
                data(11))};	
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
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the rows of parameters for the gravity loading and Save them to
% Model data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (Model.NumGravityNodes ~= 0)
    fgets(InputFile);                   % skip 2 lines
    fgets(InputFile);
    data=zeros(1,3);
    for i = 1:Model.NumGravityNodes
        data(i,:) = str2num(fgets(InputFile));            
    end
    Model = CreateStaticLoading(Model, data);
end
fgets(InputFile);                   % skip 2 lines
fgets(InputFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the parameters for the hybrid testing parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = str2num(fgets(InputFile));    
% Add the periods and damping ratio to the structure
Model.T1 = data(1);
Model.T2 = data(2);
Model.DampingRatio = data(3);
fgets(InputFile);                   % skip 2 lines
fgets(InputFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the parameters for the integration algorithm and creat integrator
% data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data2 = str2num(fgets(InputFile));
% Create a new integrator with the EQ scale, timestep, Interpolations, 
% Integration method ID and parameters
Integrator = CreateIntegrator(data(4), data(5), data(6), data2(1), data2(2));
fclose(InputFile);

% temporary code for checking model configuration visually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the model based on the information 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=0;
j=0;
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

figure (1)
line (X,Y,'Color',[0,0,1],'LineWidth',[0.5]);
hold on
% for i=1:k
%     
%     w = abs(Nodes(PZNd(i,3)).Xcoord- Nodes(PZNd(i,1)).Xcoord);
%     h = abs(Nodes(PZNd(i,4)).Ycoord- Nodes(PZNd(i,2)).Ycoord);
%     
%     PX(1,i) = Nodes(PZNd(i,1)).Xcoord;
%     PX(2,i) = Nodes(PZNd(i,2)).Xcoord+w/2;
%     PX(3,i) = Nodes(PZNd(i,3)).Xcoord;
%     PX(4,i) = Nodes(PZNd(i,4)).Xcoord-w/2;
%     PY(1,i) = Nodes(PZNd(i,1)).Ycoord-h/2;
%     PY(2,i) = Nodes(PZNd(i,2)).Ycoord;
%     PY(3,i) = Nodes(PZNd(i,3)).Ycoord+h/2;
%     PY(4,i) = Nodes(PZNd(i,4)).Ycoord;
%     PX(5,i) = PX(1,i);
%     PY(5,i) = PY(2,i);     
% end
% if num2int(k) < 1
%     line (PX,PY,'Color',[0,0,1],'LineWidth',[0.5]);
%     hold on;    
% end

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

end
  
function [par matID] = SectionInfo(data, Materials)

% search for the Material ID assigned to the section from material list
secType= data(2);
switch secType
    case 1 % wide flange section
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
