%% Update acceleration vector for KR Explicit alpha method
function Integrator = UpdateAccelerationKRExplicit(Integrator, Structure, step)

IntmPEFF = (1 - Integrator.Alphaf) * Integrator.PEFF(step,:)' + ...
                Integrator.Alphaf * Integrator.PEFF(step-1,:)';

Integrator.Acceleration = Integrator.MassMatrixHat1Inv *...
                          (Integrator.EQScaleFactor * IntmPEFF - Integrator.IntmInherentDampingForce - ...
                          Integrator.IntmRestoringForce - Integrator.PrevWeightedInertiaForce + ...
                          Structure.P0);

