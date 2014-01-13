%% Set parameters and initial conditions for the Integrator
function Integrator = InitializeIntegrator(Integrator, Structure)

% Set integration parameters depending on Integration method
switch (Integrator.MethodID)
    case 1 % CR
        % Integrator.Alpha1 = (4 * Structure.MassMatrixFree + ...
        %                     2 * Structure.DampingMatrixFreeIntegratorSetup * Integrator.Timestep + ...
        %                     Structure.StiffnessMatrixFree * Integrator.Timestep^2)^(-1) * 4 * Structure.MassMatrixFree;
        % Integrator.Alpha2 = Integrator.Alpha1;              
        alpha_1_n1 = Integrator.Lambda^2 + 2*Integrator.Lambda + 1;
        alpha_2_n1 = Integrator.Lambda + 1;
        alpha_m = 2*Integrator.Lambda^2 + 4*Integrator.Lambda + 2;
        alpha_c = 3 + 2*Integrator.Lambda - Integrator.Lambda^2;
        alpha_k = 2;
        Integrator.Alpha1 = (alpha_m * Structure.MassMatrixFree + ...
                            alpha_c * Structure.DampingMatrixFreeIntegratorSetup * Integrator.Timestep + ...
                            alpha_k * Structure.StiffnessMatrixFree * Integrator.Timestep^2)^(-1) * 2*alpha_1_n1 * Structure.MassMatrixFree;
        Integrator.Alpha2 = (alpha_m * Structure.MassMatrixFree + ...
                            alpha_c * Structure.DampingMatrixFreeIntegratorSetup * Integrator.Timestep + ...
                            alpha_k * Structure.StiffnessMatrixFree * Integrator.Timestep^2)^(-1) * 4*alpha_2_n1 * Structure.MassMatrixFree;
end

% Initial array conditions
Integrator.Displacement = zeros(Structure.NumFreeDOF, 1);
Integrator.Velocity = zeros(Structure.NumFreeDOF, 1);
%Integrator.Acceleration = Structure.MassMatrixFreeInv * Integrator.PEFF(1,:)' * Integrator.EQScaleFactor;
Integrator.Acceleration = zeros(Structure.NumFreeDOF, 1);
Integrator.RestoringForce = zeros(Structure.NumFreeDOF, 1);
