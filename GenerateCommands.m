%% Generate next integration displacement command and velocity vectors
function Integrator = GenerateCommands(Integrator, Structure, step)

% Generate the displacement and velocity vectors depending on integration methods chosen
switch (Integrator.MethodID)
    case 1 % CR
        % Displacement vector
        Integrator.Displacement = Integrator.Displacement + ...
                                  Integrator.Timestep * Integrator.Velocity + ...
                                  Integrator.Alpha2 * Integrator.Timestep^2 * Integrator.Acceleration;       
                              
        % Velocity vector
        Integrator.Velocity = Integrator.Velocity + ...
                              Integrator.Alpha1 * Integrator.Timestep * Integrator.Acceleration;
    case 2 % Rosenbrock-W
        % First Stage of RW for 1/2 the timestep
        if (Integrator.Stage == 1)  
            % Subtract Structure.P0 from the effective earthquake force          
            Integrator.e_tilda = Structure.MassMatrixFreeInv * Integrator.Timestep * ...
                                (Integrator.EQScaleFactor * Integrator.PEFF(step,:)' - ...
                                (Integrator.RestoringForce -Structure.P0 +  Structure.DampingMatrixFree * Integrator.Velocity) + ...
                                Integrator.Gamma * Integrator.Timestep * (-Structure.StiffnessMatrixFree * Integrator.Ud));
            Integrator.d_tilda = Integrator.Timestep * ...
                                (Integrator.Ud + Integrator.Gamma * Integrator.e_tilda);            
            Integrator.Displacement = Integrator.U + 1/2 * Integrator.d_tilda;  
            Integrator.Velocity = Integrator.Ud + 1/2 * Integrator.e_tilda;    
        end            
        % Second Stage of RW for 1/2 the timestep
        if (Integrator.Stage == 2)   
	    % Subtract Structure.P0 from the effective earthquake force                   
            e =  Structure.MassMatrixFreeInv * Integrator.Timestep * ...
                    (Integrator.EQScaleFactor * Integrator.PEFF(step,:)' - ...
                    (Integrator.RestoringForce -Structure.P0 + Structure.DampingMatrixFree * Integrator.Velocity) + ...
                    (Integrator.Gamma * Integrator.Alpha1 * Integrator.Timestep * Structure.StiffnessMatrixFree + ...
                    Integrator.Gamma * Structure.DampingMatrixFree) * Integrator.e_tilda);
            d =  Integrator.Timestep * ...
                    (Integrator.Ud + Integrator.Alpha2 * Integrator.e_tilda + ...
                    Integrator.Gamma * e);
            Integrator.Displacement = Integrator.U + d;              
            Integrator.Velocity = Integrator.Ud + e;              
            Integrator.U = Integrator.Displacement;
            Integrator.Ud = Integrator.Velocity;
        end   
        % Switch Stages for Rosenbrock-W method
        if (Integrator.Stage == 1)
            Integrator.Stage = 2;
        elseif (Integrator.Stage == 2)
            Integrator.Stage = 1;
        else
        end
end
                        
            
            