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
    case 2 % Rosenbrock-W        
        error('Rosenbrock-W support has been removed');
    case 3 % KR Semi-explicit alpha-method
        Alpha_Den = Structure.MassMatrixFree + Integrator.Gamma * Integrator.Timestep * Structure.DampingMatrixFreeIntegratorSetup + ...
            Integrator.Beta * Integrator.Timestep^2 * Structure.StiffnessMatrixFree;
        
        Integrator.Alpha1 = Alpha_Den^(-1) * (Structure.MassMatrixFree + ...
            Integrator.Gamma * Integrator.Timestep * Structure.DampingMatrixFreeIntegratorSetup);
        
        Integrator.Alpha2 = 0.5 * Alpha_Den^(-1) * (Structure.MassMatrixFree + ...
            (Integrator.Gamma - 2 * Integrator.Beta) * Integrator.Timestep * Structure.DampingMatrixFreeIntegratorSetup);
        
        Integrator.Alpha3 = Alpha_Den^(-1) * (Integrator.Alpham * (Structure.MassMatrixFree + ...
                            Integrator.Gamma * Integrator.Timestep * Structure.DampingMatrixFreeIntegratorSetup) + ...
                            Integrator.Alphaf * Integrator.Beta * Integrator.Timestep^2 * Structure.StiffnessMatrixFree);
        
        Integrator.MassMatrixHat1 = Structure.MassMatrixFree * (eye(size(Structure.MassMatrixFree)) - ...
            Integrator.Alpha3) + (1 - Integrator.Alphaf) * Integrator.Gamma * Integrator.Timestep * Structure.DampingMatrixFree;
        Integrator.MassMatrixHat1Inv = inv(Integrator.MassMatrixHat1);
        
        Integrator.MassMatrixHat2 = Structure.MassMatrixFree * Integrator.Alpha3 + ...
            (1 - Integrator.Alphaf) * (1 - Integrator.Gamma) * Integrator.Timestep * Structure.DampingMatrixFree;

    case 4 % KR Explicit alpha-method
        Alpha_Den = Structure.MassMatrixFree + Integrator.Gamma * Integrator.Timestep * Structure.DampingMatrixFreeIntegratorSetup + ...
            Integrator.Beta * Integrator.Timestep^2 * Structure.StiffnessMatrixFree;
        
        Integrator.Alpha1 = Alpha_Den^(-1) * Structure.MassMatrixFree;
        
        Integrator.Alpha2 = 0.5 * (1 + 2 * Integrator.Gamma) * Integrator.Alpha1;
        
        Integrator.Alpha3 = Alpha_Den^(-1) * (Integrator.Alpham * Structure.MassMatrixFree + ...
                            Integrator.Alphaf * Integrator.Gamma * Integrator.Timestep * Structure.DampingMatrixFreeIntegratorSetup + ...
                            Integrator.Alphaf * Integrator.Beta * Integrator.Timestep^2 * Structure.StiffnessMatrixFree);
        Integrator.MassMatrixHat1 = Structure.MassMatrixFree * (eye(size(Structure.MassMatrixFree)) - Integrator.Alpha3);
        Integrator.MassMatrixHat1Inv = inv(Integrator.MassMatrixHat1);
        Integrator.MassMatrixHat2 = Structure.MassMatrixFree * Integrator.Alpha3;
    case 5 % Chang 2 Int Para
        Alpha_Den = eye(size(Structure.MassMatrixFree)) + 0.5 * Integrator.Timestep * Structure.MassMatrixFreeInv *...
                    Structure.DampingMatrixFreeIntegratorSetup + 0.25 * Integrator.Timestep^2 * Structure.MassMatrixFreeInv * Structure.StiffnessMatrixFree;
                
        Integrator.Alpha1 = Alpha_Den^(-1) * (eye(size(Structure.MassMatrixFree)) + ...
            0.5 * Integrator.Timestep * Structure.MassMatrixFreeInv * Structure.DampingMatrixFreeIntegratorSetup);
        
        Integrator.Alpha2 = 0.5 * Alpha_Den^(-1);
    case 6 % Chang 3 Int Para
        Alpha_Den = eye(size(Structure.MassMatrixFree)) + 0.5 * Integrator.Timestep * Structure.MassMatrixFreeInv *...
                    Structure.DampingMatrixFreeIntegratorSetup + 0.25 * Integrator.Timestep^2 * Structure.MassMatrixFreeInv * Structure.StiffnessMatrixFree;
                
        Integrator.Alpha1 = Alpha_Den^(-1) * (eye(size(Structure.MassMatrixFree)) + ...
            0.5 * Integrator.Timestep * Structure.MassMatrixFreeInv * Structure.DampingMatrixFreeIntegratorSetup);
        
        Integrator.Alpha2 = Alpha_Den^(-1) * (eye(size(Structure.MassMatrixFree)) + ...
            0.25 * Integrator.Timestep * Structure.MassMatrixFreeInv * Structure.DampingMatrixFreeIntegratorSetup);
        
        Integrator.Alpha3 = 0.25 * Alpha_Den^(-1);
end

% Initial array conditions
Integrator.Displacement = zeros(Structure.NumFreeDOF, 1);
Integrator.Velocity = zeros(Structure.NumFreeDOF, 1);
Integrator.Acceleration = Structure.MassMatrixFreeInv * Integrator.PEFF(1,:)' * Integrator.EQScaleFactor;
Integrator.RestoringForce = zeros(Structure.NumFreeDOF, 1);
