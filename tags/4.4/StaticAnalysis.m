%% Calculate matrices and return back the Structure object with new matricies
function [Structure, Elements, Integrator, flag] = StaticAnalysis(Structure, Elements, Integrator)

% Analysis 
Integrator.Displacement = Structure.StiffnessMatrixFree \ Structure.P0;

for element = 1:Structure.NumElements
    if Elements{element}.Type == 9
        ndof = 12;
        % Get the resisting force of each element.
        [Elements{element}, ElementRF] = ...
            Elements{element}.GetRestoringForce(Elements{element}, Structure, Integrator);
        % Assemble the structure restoring force vector based on the
        % degrees of freedom node locations.
        dof=zeros(ndof,1);
        dof(1,1) =Elements{element}.Nodes(1).UX;
        dof(2,1) =Elements{element}.Nodes(1).UY;
        dof(3,1) =Elements{element}.Nodes(1).THETA;
        dof(4,1) =Elements{element}.Nodes(2).UX;
        dof(5,1) =Elements{element}.Nodes(2).UY;
        dof(6,1) =Elements{element}.Nodes(2).THETA;
        dof(7,1) =Elements{element}.Nodes(3).UX;
        dof(8,1) =Elements{element}.Nodes(3).UY;
        dof(9,1) =Elements{element}.Nodes(3).THETA;
        dof(10,1)=Elements{element}.Nodes(4).UX;
        dof(11,1)=Elements{element}.Nodes(4).UY;
        dof(12,1)=Elements{element}.Nodes(4).THETA;
        
    else
        ndof = 6;
        % Get the resisting force of each element.
        [Elements{element}, ElementRF] = ...
            Elements{element}.GetRestoringForce(Elements{element}, Structure, Integrator);
        % Assemble the structure restoring force vector based on the
        % degrees of freedom node locations.
        dof=zeros(ndof,1);
        dof(1,1)=Elements{element}.Nodes(1).UX;
        dof(2,1)=Elements{element}.Nodes(1).UY;
        dof(3,1)=Elements{element}.Nodes(1).THETA;
        dof(4,1)=Elements{element}.Nodes(2).UX;
        dof(5,1)=Elements{element}.Nodes(2).UY;
        dof(6,1)=Elements{element}.Nodes(2).THETA;
    end
    
    for i=1:ndof
        if dof(i,1) ~= -1
            Integrator.RestoringForce(dof(i,1), 1) = ElementRF(i, 1) + Integrator.RestoringForce(dof(i,1), 1);
        end
    end
end

% check equilibrium
resP = Integrator.RestoringForce - Structure.P0;
nresidual = norm(resP);

if nresidual < 0.0000001
    flag =1;
else flag = 0;
end

