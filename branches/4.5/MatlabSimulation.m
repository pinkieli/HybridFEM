%% Dynamic steps of the integration algorithm for the Matlab Simulation
for step = 2:Integrator.Steps
    Am_Out_Steps(step) = step;
    fprintf('Step %i/%i\n',step,Integrator.Steps);
    % Generate Commands and store them in the Results variable.
    Integrator = GenerateCommands(Integrator, Structure, step);
    Am_Out_Displacement(step,:) = Integrator.Displacement';
    switch (Integrator.MethodID)
        case 1 % CR
            Am_Out_Velocity(step,:)=Integrator.Velocity';
        case 2 % Rosenbrock-W
            Am_Out_Velocity(step,:)=Integrator.Velocity';
        case 3 % KR Semi-explicit alpha-method
            EffForcePrevStep = Integrator.MassMatrixHat2 * Integrator.Acceleration + ...
                Structure.DampingMatrixFree * Integrator.Velocity + Integrator.Alphaf * ...
                (Integrator.RestoringForce - Integrator.EQScaleFactor * Integrator.PEFF((step - 1),:)');
        case 4 % KR Explicit alpha-method
            Am_Out_Velocity(step,:)=Integrator.Velocity';
            Integrator.IntmInherentDampingForce = Structure.DampingMatrixFree * ...
                ((1 - Integrator.Alphaf) * Am_Out_Velocity(step,:)' + ...
                Integrator.Alphaf * Am_Out_Velocity(step-1,:)');
            Integrator.PrevWeightedInertiaForce = Integrator.MassMatrixHat2 * Integrator.Acceleration;
    end
    % Calculate Restoring Forces per Element and form full structure
    % restoring force vector.
    Integrator.RestoringForce = zeros(Structure.NumFreeDOF, 1);   
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
        
        % Record element restoring forces
%         if exist('ElementRestoringForce','var') == 1
         if isnumeric(Structure.RecordElements) == 1
            for i = 1 : length(Structure.RecordElements)
                if Elements{element}.ID == Am_Out_ElementRestoringForce(1,i).ID;
                    Am_Out_ElementRestoringForce(1,i).Resp(step,:) = ElementRF;
                end
            end
        end
    end
    % Store restoring force results.
    Am_Out_RestoringForce(step,:) = Integrator.RestoringForce';    
    
    % Update Acceleration Vector
    switch (Integrator.MethodID)
        case 1 % CR
            Integrator = UpdateAccelerationCR(Integrator, Structure, step);
            Am_Out_Acceleration(step,:)=Integrator.Acceleration';
            
        case 3 % KR Semi-explicit alpha-method            
            Integrator = UpdateAccelerationKRSemiExplicit(Integrator, Structure, step, EffForcePrevStep);
            Am_Out_Acceleration(step,:)=Integrator.Acceleration';
            Integrator = UpdateVelocityKRSemiExplicit(Integrator, Am_Out_Acceleration(step-1:step,:));
            Am_Out_Velocity(step,:)=Integrator.Velocity';
            
        case 4 % KR Explicit alpha-method
            Integrator.IntmRestoringForce = (1 - Integrator.Alphaf) * Am_Out_RestoringForce(step,:)' + ...
                Integrator.Alphaf * Am_Out_RestoringForce(step-1,:)';
            Integrator = UpdateAccelerationKRExplicit(Integrator, Structure, step);
            Am_Out_Acceleration(step,:)=Integrator.Acceleration';
            
        case 5 % Chang 2 Int Para
            Integrator = UpdateVelocityChang(Integrator, Structure, step);
            Am_Out_Velocity(step,:)=Integrator.Velocity';
            Integrator = UpdateAccelerationCR(Integrator, Structure, step);
            Am_Out_Acceleration(step,:)=Integrator.Acceleration';
            
        case 6 % Chang 3 Int Para
            Integrator = UpdateVelocityChang(Integrator, Structure, step);
            Am_Out_Velocity(step,:)=Integrator.Velocity';
            Integrator = UpdateAccelerationCR(Integrator, Structure, step);
            Am_Out_Acceleration(step,:)=Integrator.Acceleration';
    end
    %figure(1)
    %dof=1; 
    %plot(DisplacementResults([1:step],dof),RestoringForceResults([1:step],dof),'r',DisplacementResults(step,dof),RestoringForceResults(step,dof),'--b>','MarkerSize',5,'MarkerFaceColor','b','LineWidth',2);
    %figure(2)
    %STEP=[1:step]';
    %plot(STEP([1:step],1),DisplacementResults([1:step],dof),'r',STEP(step,1),DisplacementResults(step,dof),'--b>','MarkerSize',5,'MarkerFaceColor','b','LineWidth',2);
    %figure(3)
    %plot(STEP([1:step],1),RestoringForceResults([1:step],dof),'r',STEP(step,1),RestoringForceResults(step,dof),'--b>','MarkerSize',5,'MarkerFaceColor','b','LineWidth',2);
pause(0.0005);
end

%% Plot global results
time = (0:1:Integrator.Steps-1)'*Integrator.Timestep;

figure;
plot(time,Am_Out_Velocity);
title('Velocity');

figure;
plot(time,Am_Out_RestoringForce);
title('Force');
 
figure;
plot(time,Am_Out_Acceleration);
title('Accelerations');

figure; 
plot(time, Am_Out_Displacement);
title('Displacement');