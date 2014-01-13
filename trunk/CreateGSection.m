%% Create General Section object (not used yet, planned to use with new input reader)
function section = CreateGSection(id, type, data, mat, UICheck)
section.ID = id;
section.Type = type;
if type ==3
    Nf = length(data(:,1));
    
    fiber.yi   = [ 0 0 ];
    fiber.A    =  0;
    fiber.prop = mat{1,data(1,3)};
    TotalArea=0;
    
    for i=1:Nf
        yf=data(i,1);
        fiber.yi   = [ yf 1];
        fiber.A    =  data(i,2);
        TotalArea=TotalArea+fiber.A;
        fiber.prop = mat{1,data(i,3)};
        fibers(i) = fiber;
    end
    
    if (UICheck == 1)
        fprintf('[FEM] Check point: Total area of the general section is: %f\n', TotalArea);    
        choice = questdlg('Is the total area of the general section correct?','Section Checking', ...
            'No', 'Yes','Yes');
        if strcmp(choice,'No')
            fprintf('Return and check your general fiber section');
        end
    end
else
    msg = ['Unknown Section Type ' str2num(section.Type)];
    errordlg(msg,'Input Error');
end

% Create section object
section.fibers = fibers;	% Fiber data
section.numfibers = length(fibers);
section.vs	   = [0 0]';	% Section deformations
