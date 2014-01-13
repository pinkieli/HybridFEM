%% Create an element forces recorder 
function ElementRestoringForce = CreateElementRecorder(Structure, Elements, Integrator)
NumElements = length(Structure.RecordElements);
if (NumElements > 0)
    for i = NumElements : -1 : 1 %for implicit preallocation
        if Elements{Structure.RecordElements(i)}.Type == 9
            ndof = 12;
        else
            ndof = 6;
        end
        Out_ElementRestoringForce(i).ID = Structure.RecordElements(i);
        Out_ElementRestoringForce(i).Resp = zeros(Integrator.Steps,ndof);    
    end    
    ElementRestoringForce = Out_ElementRestoringForce;
else
    ElementRestoringForce = 0;
end


