%% Update velocity vector for KR Semi-explicit alpha method
function Integrator = UpdateVelocityKRSemiExplicit(Integrator,Acceleration)

Integrator.Velocity = Integrator.Velocity + Integrator.Timestep * ((1 - Integrator.Gamma) * Acceleration(1,:)' + ...
                      Integrator.Gamma * Acceleration(2,:)');
