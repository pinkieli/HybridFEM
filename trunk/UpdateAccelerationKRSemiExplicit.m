%% Update acceleration vector for KR Semi-explicit alpha method
function Integrator = UpdateAccelerationKRSemiExplicit(Integrator, Structure, step, EffForcePrevStep)

% Subtract Structure.P0 attributed to gravity load from effective earthquake force
Integrator.Acceleration = Integrator.MassMatrixHat1Inv * ...
                          ((1 - Integrator.Alphaf) * (Integrator.EQScaleFactor * Integrator.PEFF(step,:)' - ...
                          Integrator.RestoringForce) - EffForcePrevStep + Structure.P0);
                    