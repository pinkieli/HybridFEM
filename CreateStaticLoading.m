%% Calculate matrices and return back the Structure object with new matricies
function Structure = CreateStaticLoading(Structure, data)

% initial condition for the time history analysis
Structure.P0 = zeros(Structure.NumFreeDOF,1);
for i=1:size(data,1)
    dofNum(i) = Structure.DOF(data(i,2),data(i,1));
end
Structure.P0(dofNum) = data(:,3);

