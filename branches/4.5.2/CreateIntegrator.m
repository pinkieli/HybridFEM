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
        error('Rosenbrock-W support has been removed');
    case 3 % KR Semi-explicit alpha-method
        Rhoinfy = Param1;
        vars.Alpham = (2 * Rhoinfy - 1) / (Rhoinfy + 1);
        vars.Alphaf = Rhoinfy / (Rhoinfy + 1);
        vars.Gamma = 0.5 - vars.Alpham + vars.Alphaf;
        vars.Beta = 0.25 * (0.5 + vars.Gamma)^2;
    case 4 % KR Explicit alpha-method
        Rhoinfy = Param1;
        vars.Alpham = (2 * Rhoinfy - 1) / (Rhoinfy + 1);
        vars.Alphaf = Rhoinfy / (Rhoinfy + 1);
        vars.Gamma = 0.5 - vars.Alpham + vars.Alphaf;
        vars.Beta = 0.25 * (0.5 + vars.Gamma)^2;
    case 5 % Chang 2 Int Para
        vars.Gamma = 0;  % placeholder for Simulink compatibility        
        vars.Lambda = Param1; % Not used
    case 6 % Chang 3 Int Para
        vars.Gamma = 0;  % placeholder for Simulink compatibility        
        vars.Lambda = Param1; % Not used
    otherwise 
        error('No Integration method selected');
end