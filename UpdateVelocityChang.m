%% Update velocity vector for Chang's algorithms
function Integrator = UpdateVelocityChang(Integrator, Structure, step)

Integrator.Velocity = Structure.MassDampingInv * (Structure.MassMatrixFree * (Integrator.Velocity + 0.5 * Integrator.Timestep * Integrator.Acceleration) + ...
                      0.5 * Integrator.Timestep * (Integrator.EQScaleFactor * Integrator.PEFF(step,:)' - Integrator.RestoringForce + Structure.P0));
