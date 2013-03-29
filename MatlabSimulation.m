%% Dynamic steps of the integration algorithm for the Matlab Simulation
figure;
for step = 1:Integrator.Steps    
    Out_Steps(step) = step;
%     if (mod(step,100) == 0)
        fprintf('Step %i/%i\n',step,Integrator.Steps);
%     end
    % Generate Commands and store them in the Results variable.
    Integrator = GenerateCommands(Integrator, Structure, step);    
    Out_Displacement(step,:) = Integrator.Displacement';
    Out_Velocity(step,:) = Integrator.Velocity';
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
    end
    % Store restoring force results.
    Out_RestoringForce(step,:) = Integrator.RestoringForce';  
    
    % Update Acceleration Vector only for CR Method
    if (Integrator.MethodID == 1)
        Integrator = UpdateAccelerationCR(Integrator, Structure, step);
        Out_Acceleration(step,:)=Integrator.Acceleration';    	 
    end            
end

%% Plot global results
time = (0:1:Integrator.Steps-1)'*Integrator.Timestep;

figure(2)
plot(time,Out_Velocity);
title('Velocity');

figure(3)
plot(time,Out_RestoringForce);
title('Restoring Force');
 
figure(4)
plot(time,Out_Acceleration);
title('Acceleration');

figure(5); 
plot(time, Out_Displacement);
title('Displacement');