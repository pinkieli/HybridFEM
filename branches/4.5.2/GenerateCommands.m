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
    case 3 % KR Semi-explicit alpha-method
        % Displacement vector
        Integrator.Displacement = Integrator.Displacement + ...
                                  Integrator.Alpha1 * Integrator.Timestep * Integrator.Velocity + ...
                                  Integrator.Alpha2 * Integrator.Timestep^2 * Integrator.Acceleration;
                              
    case 4 % KR Explicit alpha-method
        % Displacement vector
        Integrator.Displacement = Integrator.Displacement + ...
                                  Integrator.Timestep * Integrator.Velocity + ...
                                  Integrator.Alpha2 * Integrator.Timestep^2 * Integrator.Acceleration;
        % Velocity vector
        Integrator.Velocity = Integrator.Velocity + ...
                              Integrator.Alpha1 * Integrator.Timestep * Integrator.Acceleration;
                          
    case 5 % Chang 2 Int Para
        % Displacement vector
        Integrator.Displacement = Integrator.Displacement + ...
                                  Integrator.Alpha1 * Integrator.Timestep * Integrator.Velocity + ...
                                  Integrator.Alpha2 * Integrator.Timestep^2 * Integrator.Acceleration;
                              
    case 6 % Chang 3 Int Para
        % Displacement vector
        Integrator.Displacement = Integrator.Alpha1 * Integrator.Displacement + ...
                                  Integrator.Alpha2 * Integrator.Timestep * Integrator.Velocity + ...
                                  Integrator.Alpha3 * Integrator.Timestep^2 * Integrator.Acceleration;
end
                        
            
            