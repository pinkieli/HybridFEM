%% Get the ground acceleration record and prepare as an effective force record
function Integrator = GetEQHistory(EQ_FILE, Integrator, Structure, Elements)

% Modify EQ History depending on Integration Method
switch (Integrator.MethodID)
    case 1 % CR
        dt = Integrator.Timestep;
    case 2 % Rosenbrock-W
        dt = Integrator.Timestep/2;
    case 3 % KR-II(b) alpha-method
        dt = Integrator.Timestep;
    case 4 % Chang 2 Int Para
        dt = Integrator.Timestep;
    case 5 % Chang 3 Int Para
        dt = Integrator.Timestep;
    case 6 % KR-II(c) alpha-method
        dt = Integrator.Timestep;
end

% Load file
EQ = load(EQ_FILE);
% Interpolate the ground motion depending on the sample rate 
gtime = EQ(:,1); % Time
xgdd = EQ(:,2);  % Record
time = 0:dt:gtime(end); time=time';
EQ = interp1(gtime,xgdd,time);

% Form influcence vector
InfluenceVector = zeros(Structure.NumFreeDOF, 1);
for i = 1:Structure.NumElements
    if (Elements{i}.Type == 9)
        for j = 1:4
            if (Elements{i}.Nodes(j).UX ~= -1)
                InfluenceVector(Elements{i}.Nodes(j).UX, 1) = 1.0;
            end
        end
    else
        for j = 1:2
            if (Elements{i}.Nodes(j).UX ~= -1)
                InfluenceVector(Elements{i}.Nodes(j).UX, 1) = 1.0;
            end
        end
    end
end

% Form Effective Force input.  No Gravity included.  
% Scaling Factor included inside Integrator.
NegativeEFF = -1.0 * Structure.MassMatrixFree * InfluenceVector;

% Save number of steps and effective force array to Integrator
Integrator.Steps = length(EQ);
Integrator.PEFF = [ EQ * NegativeEFF' ];
