%% Update acceleration vector for CR method
function Integrator = UpdateAccelerationCR(Integrator, Structure, step)

% Subtract Structure.P0 attributed to gravity load from effective earthquake force
Integrator.Acceleration = Structure.MassMatrixFreeInv * ...
                          (Integrator.EQScaleFactor * Integrator.PEFF(step,:)' - Structure.DampingMatrixFree * ...
                           Integrator.Velocity - Integrator.RestoringForce + Structure.P0);
