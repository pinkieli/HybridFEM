%% Number the DOFs
function [Model Nds]= DOFNumberer(Model, nsDOF, Nds)
% numbering degrees of freedom including constraints and restraint 
% Model = Model data structure 
% nsDOF =  number of slaved dofs
% Nds = Node data structure 

nNodes = Model.NumNodes;                % total number of nodes
BOUND  = Model.BOUND;                   % array with boundary conditions
ndfpn  = Model.NumDOF;                  % number of dof per node
EQDOF  = Model.EqDOF;                   % array with constraint conditions 

ntDOF = ndfpn*nNodes;                   % number of total dof's
nrDOF = nnz(BOUND);                     % number of nonzero elements in BOUND = number of restrained dof's
nfDOF = ntDOF-nrDOF-nsDOF;              % number of free dof's  

% 
dof  = zeros(ndfpn,nNodes);             % initialize array for dof numbers
rInd = find(BOUND');                    % locate non-zero elements (restrained dof) in BOUND 
dof(rInd) = -1;                         % for nrDOFs, assign number -1 to index locations in DOF 
rsInd = rInd;                           % initialize array for restrained and constrained dofs

for k=1:max(max(EQDOF))
    eqdInd = find(EQDOF'==k);           % locate constrained dofs in EQDOF 
    mInds(k) = eqdInd(1);               % assign master node dof index 
    sInds{k} = eqdInd(2:length(eqdInd));% assign slave node dofs' indices 
    rsInd = [rsInd; sInds{k}];          % assign dof indices for restrained and constrained dof   
end

fInd = setdiff(1:ndfpn*nNodes, rsInd);  % free dof's are the difference between the total and restrained and constrained dofs indices
dof(fInd) = 1:nfDOF;                    % assign numbers 1 to nfDOF to indx locations in array dof for free dof's

% fill the slave dof number (zero elements in dof) with the master node dof number 
for k=1:max(max(EQDOF))
    dof(sInds{k}) =dof(mInds(k));     
end

% check model data structure
if ~(Model.NumFreeDOF == nfDOF)
    msg =['Check the number of free, restrained and slaved dofs'];
    errordlg(msg,'Input Error');
end

% update Model dof numbers
Model.DOF  = dof;

% assign dof numbers including constraints to node data structure
for i=1:Model.NumNodes
    Nds(i).UX =Model.DOF(1,i);
    Nds(i).UY =Model.DOF(2,i);
    Nds(i).THETA =Model.DOF(3,i);
end

