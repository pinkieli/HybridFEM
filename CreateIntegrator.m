%% Create an Integrator
function vars = CreateIntegrator(EQScaleFactor, Timestep, Interpolations, MethodID, Param1)

% Set integration parameters
vars.EQScaleFactor = EQScaleFactor;
vars.Timestep = Timestep;
vars.TimestepSquared = Timestep^2;
vars.MethodID = MethodID;
vars.Interpolations = Interpolations;

% Definable parameters for the selected integration method
switch (MethodID)
    case 1 % CR
        vars.Gamma = 0;  % placeholder for Simulink compatibility        
        vars.Lambda = Param1;
    case 2 % Rosenbrock-W
        vars.Gamma = Param1;        
        % Check that Interpolations are even since they need to be halved
        % for both stages
        if (mod(Interpolations,2))
            error('Interpolations need to be even for Rosenbrock-W Integrator');
        end
    otherwise 
        error('No Integration method selected');
end