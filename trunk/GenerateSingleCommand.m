%% Generate next integration displacement command and velocity vectors
function Integrator = GenerateSingleCommand(DOF, Integrator, Structure, step)

    
if step == 2
    Integrator.SingleCommand = 0;
end

        % Displacement vector        
        Integrator.SingleCommand = Integrator.SingleCommand + ...
                                  Integrator.Timestep * Integrator.Velocity(DOF) + ...
                                  Integrator.Alpha2(DOF,:) * Integrator.Timestep^2 * Integrator.Acceleration;

                              