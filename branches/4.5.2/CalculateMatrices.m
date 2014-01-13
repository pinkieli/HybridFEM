%% Calculate matrices and return back the Structure object with new matrices
% DO NOT MODIFY THIS FILE! Instead, copy it to your working folder and make
% changes there. That file will run instead of this one
function Structure = CalculateMatrices(Structure, Elements, Integrator)
      
% Preset arrays
MassMatrixFree = zeros(Structure.NumFreeDOF, Structure.NumFreeDOF);
StiffnessMatrixFree=zeros(Structure.NumFreeDOF, Structure.NumFreeDOF);
DampingMatrixFree=zeros(Structure.NumFreeDOF, Structure.NumFreeDOF);
DampingMatrixFreeIntegratorSetup=zeros(Structure.NumFreeDOF, Structure.NumFreeDOF);

% Initial damping calculations
PI = 3.141592654;
omega1 = 2.0 * PI / Structure.T1;
omega2 = 2.0 * PI / Structure.T2;
A0 = Structure.DampingRatio * 2.0 * omega1 * omega2 / (omega1 + omega2);
A1 = Structure.DampingRatio * 2.0 / (omega1 + omega2);

% Go through elements and set the degrees of freedom
for i=1:Structure.NumElements
    [ElementStiffnessMatrix,...
    ElementMassMatrix] = Elements{i}.FormElementMatrices(Elements{i}, 0);

    if (Elements{i}.Type==9) 
        ndof = 12;
        dof=zeros(ndof,1);
        dof(1,1)=Elements{i}.Nodes(1).UX;
        dof(2,1)=Elements{i}.Nodes(1).UY;
        dof(3,1)=Elements{i}.Nodes(1).THETA;
        dof(4,1)=Elements{i}.Nodes(2).UX;
        dof(5,1)=Elements{i}.Nodes(2).UY;
        dof(6,1)=Elements{i}.Nodes(2).THETA;
        dof(7,1)=Elements{i}.Nodes(3).UX;
        dof(8,1)=Elements{i}.Nodes(3).UY;
        dof(9,1)=Elements{i}.Nodes(3).THETA;
        dof(10,1)=Elements{i}.Nodes(4).UX;
        dof(11,1)=Elements{i}.Nodes(4).UY;
        dof(12,1)=Elements{i}.Nodes(4).THETA; 
        
        if( Structure.RigidLinkNo )
            TR = eye(12,12); % Transformation matrix for rigid link            
            for jjj=1:4
                idx=find(Structure.RigidLinkNodeID(:,2) == Elements{i}.Nodes(jjj).ID); % find if a slave node exists in the element
                if ( idx )
                    TL=CalculateRigidLinkTransformation(Structure.RigidLinkMaster{idx},Structure.RigidLinkSlave{idx} );
                    TR(3*(jjj-1)+1:3*jjj,3*(jjj-1)+1:3*jjj) = TL;
                end
            end               
            ElementStiffnessMatrix = TR' * ElementStiffnessMatrix * TR;
        end
        
    else       
        ndof = 6;
        dof=zeros(ndof,1);
        dof(1,1)=Elements{i}.Nodes(1).UX;
        dof(2,1)=Elements{i}.Nodes(1).UY;
        dof(3,1)=Elements{i}.Nodes(1).THETA;
        dof(4,1)=Elements{i}.Nodes(2).UX;
        dof(5,1)=Elements{i}.Nodes(2).UY;
        dof(6,1)=Elements{i}.Nodes(2).THETA;
        
        if( Structure.RigidLinkNo )
            TR = eye(6,6); % Transformation matrix for rigid link            
            for jjj=1:2
                idx=find(Structure.RigidLinkNodeID(:,2) == Elements{i}.Nodes(jjj).ID); % find if a slave node exists in the element
                if ( idx )
                    TL=CalculateRigidLinkTransformation(Structure.RigidLinkMaster{idx},Structure.RigidLinkSlave{idx} );
                    TR(3*(jjj-1)+1:3*jjj,3*(jjj-1)+1:3*jjj) = TL;
                end
            end               
            ElementStiffnessMatrix = TR' * ElementStiffnessMatrix * TR;
        end
        
        % Experimental element Type 2 setup
        if (Elements{i}.Type==2) 
           ElementDampingMatrix=zeros(ndof,ndof);
           ElementDampingMatrix(1,1)=Elements{i}.Cexp;
           ElementDampingMatrix(1,4)=-Elements{i}.Cexp;
           ElementDampingMatrix(4,4)=Elements{i}.Cexp;
           ElementDampingMatrix(4,1)=-Elements{i}.Cexp;
        end
    end
    
    % Go through each dof and set the mass and stiffness matrices
    for ii=1:ndof
        for jj=1:ndof
            if (dof(ii,1) ~= -1  && dof(jj,1) ~= -1)            
                StiffnessMatrixFree(dof(ii,1),dof(jj,1))=StiffnessMatrixFree(dof(ii,1),dof(jj,1))+...
                                                            ElementStiffnessMatrix(ii,jj);
                MassMatrixFree(dof(ii,1),dof(jj,1))=MassMatrixFree(dof(ii,1),dof(jj,1))+...
                                                        ElementMassMatrix(ii,jj);
                if (~(Elements{i}.Type==2)) % element type 2 does not contibute to DampingMatrixFree                
                    DampingMatrixFree(dof(ii,1),dof(jj,1))=DampingMatrixFree(dof(ii,1),dof(jj,1))+...
                            Elements{i}.DampMassFac*A0*ElementMassMatrix(ii,jj)+...
                            Elements{i}.DampStiffFac*A1*ElementStiffnessMatrix(ii,jj);
                end
                DampingMatrixFreeIntegratorSetup(dof(ii,1),dof(jj,1))=DampingMatrixFree(dof(ii,1),dof(jj,1));
                if(Elements{i}.Type==2) % element type 2 contributes to DampingMatrixFreeIntegratorSetup
                    DampingMatrixFreeIntegratorSetup(dof(ii,1),dof(jj,1))=DampingMatrixFreeIntegratorSetup(dof(ii,1),dof(jj,1))+...
                                                                            ElementDampingMatrix(ii,jj);
                end
            end
        end
    end
end

%% Modified by Akbar for Experimental Element Damping contribution
DampingMatrixFreeIntegratorSetup=DampingMatrixFree;
for i=1:Structure.NumElements
    ndof = 6;
    dof=zeros(ndof,1);
    dof(1,1)=Elements{i}.Nodes(1).UX;
    dof(2,1)=Elements{i}.Nodes(1).UY;
    dof(3,1)=Elements{i}.Nodes(1).THETA;
    dof(4,1)=Elements{i}.Nodes(2).UX;
    dof(5,1)=Elements{i}.Nodes(2).UY;
    dof(6,1)=Elements{i}.Nodes(2).THETA;
    if (Elements{i}.Type==2) 
        ElementDampingMatrix=zeros(ndof,ndof);
        ElementDampingMatrix(1,1)=Elements{i}.Cexp;
        ElementDampingMatrix(1,4)=-Elements{i}.Cexp;
        ElementDampingMatrix(4,4)=Elements{i}.Cexp;
        ElementDampingMatrix(4,1)=-Elements{i}.Cexp;
        Elements{i}.Cexp		   
    end
        
    for ii=1:ndof
        for jj=1:ndof
            if (dof(ii,1) ~= -1  && dof(jj,1) ~= -1)   
                if(Elements{i}.Type==2) % element type 2 contributes to DampingMatrixFreeIntegratorSetup
				    DampingMatrixFreeIntegratorSetup(dof(ii,1),dof(jj,1))
                    DampingMatrixFreeIntegratorSetup(dof(ii,1),dof(jj,1))=DampingMatrixFreeIntegratorSetup(dof(ii,1),dof(jj,1))+...
                                                                            ElementDampingMatrix(ii,jj);			    				
                end
            end
        end
    end
end
%% Add experimental Rayleigh damping matrix (if any) to numerical Rayleigh damping matrix and damping matrix for integration parameters
if isfield(Structure,'CexpRayleigh')
    for i = 1 : length(Structure.InterfaceDOF)
        for j = 1 : length(Structure.InterfaceDOF)
            DampingMatrixFree(Structure.InterfaceDOF(i),Structure.InterfaceDOF(j)) = DampingMatrixFree(Structure.InterfaceDOF(i),Structure.InterfaceDOF(j)) + Structure.CexpRayleigh(i,j);
            DampingMatrixFreeIntegratorSetup(Structure.InterfaceDOF(i),Structure.InterfaceDOF(j)) = DampingMatrixFreeIntegratorSetup(Structure.InterfaceDOF(i),Structure.InterfaceDOF(j)) + Structure.CexpRayleigh(i,j);
        end
    end
end
%% Add experimental Initial damping matrix(if any) to the damping matrix for integration parameters
if isfield(Structure,'CexpInitial')
    for i = 1 : length(Structure.InterfaceDOF)
        for j = 1 : length(Structure.InterfaceDOF)
            DampingMatrixFreeIntegratorSetup(Structure.InterfaceDOF(i),Structure.InterfaceDOF(j)) = DampingMatrixFreeIntegratorSetup(Structure.InterfaceDOF(i),Structure.InterfaceDOF(j)) + Structure.CexpInitial(i,j);
        end
    end
end
%% Add experimental stiffness matrix (if any) to the total stiffness matrix
if isfield(Structure,'KexpInitial')
    for i = 1 : length(Structure.InterfaceDOF)
        for j = 1 : length(Structure.InterfaceDOF)
            StiffnessMatrixFree(Structure.InterfaceDOF(i),Structure.InterfaceDOF(j)) = StiffnessMatrixFree(Structure.InterfaceDOF(i),Structure.InterfaceDOF(j)) + Structure.KexpInitial(i,j);
        end
    end
end
%% Save matrices back to the structure
Structure.StiffnessMatrixFree = StiffnessMatrixFree;
Structure.MassMatrixFree = MassMatrixFree;
Structure.DampingMatrixFree = DampingMatrixFree;
Structure.DampingMatrixFreeIntegratorSetup=DampingMatrixFreeIntegratorSetup;
Structure.MassMatrixFreeInv = inv(MassMatrixFree);

%% Calculate inverse of sum of mass and damping matrices for Chang's
% algorithms
if (Integrator.MethodID == 5 || Integrator.MethodID == 6)
    Structure.MassDampingInv = inv(Structure.MassMatrixFree + 0.5 * Integrator.Timestep * Structure.DampingMatrixFree);
end