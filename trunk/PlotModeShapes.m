%% Plot mode shapes: developed by CKolay
function [] = PlotModeShapes(Structure, Elements, Nodes, NumModeShapes)
% Mode shapes are plotted only considering horizontal and vertical
% translation of nodes without considering rotations
% For higher modes user may want to change the scale multiplier below
% This factor times the corresponding dimension of the model is the max
% displacement of a node
dimMultiply = 0.15; 
%% Calculate eigenvalues and eigenvectors and sort them
[V,D] = eig(Structure.StiffnessMatrixFree,Structure.MassMatrixFree);
Dvec = zeros(length(D),1);
for i = 1 : length(D)
    Dvec(i) = D(i,i);
end
[~, IX] = sort(Dvec,'ascend');

%% Disregard Element Type 4 and 9 for plotting mode shapes
j=0;
for i=1:Structure.NumElements
    if Elements{i}.Type == 9 || Elements{i}.Type == 4
    else
        j=j+1;
        iN= Elements{i}.Nodes(1).ID;
        jN= Elements{i}.Nodes(2).ID;    
        CON{j,1} = [iN jN];
    end
end
ElmNd  = cat(1,CON{1:j,1});
%% For plotting undeformed geometry of structure
X = zeros(2, length(ElmNd)); Y = zeros(2, length(ElmNd));
NdX = zeros(1, Structure.NumNodes); NdY = zeros(1, Structure.NumNodes);
for j = 1 : length(ElmNd)
    X(1,j) = Nodes(ElmNd(j,1)).Xcoord;
    X(2,j) = Nodes(ElmNd(j,2)).Xcoord;
    Y(1,j) = Nodes(ElmNd(j,1)).Ycoord;
    Y(2,j) = Nodes(ElmNd(j,2)).Ycoord;
end
for k = 1:Structure.NumNodes
    NdX(1,k) = Nodes(k).Xcoord;
    NdY(1,k) = Nodes(k).Ycoord;
end
maxXDim = max(max(NdX)-min(NdX)); % max x-dimension of model
maxYDim = max(max(NdY)-min(NdY)); % max y-dimension of model
%% For each mode shape
MarkerSize = 25;
for i = 1 : NumModeShapes
    % Determine scale for x-axis
    Xscale = dimMultiply * maxXDim / max(abs(V(:, IX(i))));
    % Determine scale for y-axis
    Yscale = dimMultiply * maxYDim / max(abs(V(:, IX(i))));
    % for line structure with one dimension being zero
    if Xscale == 0 
        Xscale = Yscale;
    elseif Yscale == 0
        Yscale = Xscale;
    end
    figure(i+1); hold on
    line (X,Y,'Color',[0,0,1],'LineWidth',[0.5],'LineStyle','--');
    scatter (NdX,NdY,MarkerSize,'s');
    Xdef = zeros(1,length(ElmNd)); Ydef = zeros(1,length(ElmNd)); % initialize for deformed geometry
    for j = 1 : length(ElmNd)
        if Structure.DOF(1, Nodes(ElmNd(j,1)).ID) ==-1
            dX1 = 0;
        else
            dX1 = Xscale * V(Structure.DOF(1, Nodes(ElmNd(j,1)).ID), IX(i));
        end
        Xdef(1,j) = Nodes(ElmNd(j,1)).Xcoord + dX1;
        
        if Structure.DOF(1, Nodes(ElmNd(j,2)).ID) ==-1
            dX2 = 0;
        else
            dX2 = Xscale * V(Structure.DOF(1, Nodes(ElmNd(j,2)).ID), IX(i));
        end
        Xdef(2,j) = Nodes(ElmNd(j,2)).Xcoord + dX2;
        
        if Structure.DOF(2, Nodes(ElmNd(j,1)).ID) ==-1
            dY1 = 0;
        else
            dY1 = Yscale * V(Structure.DOF(2, Nodes(ElmNd(j,1)).ID), IX(i));
        end
        Ydef(1,j) = Nodes(ElmNd(j,1)).Ycoord + dY1;
        if Structure.DOF(2, Nodes(ElmNd(j,2)).ID) ==-1
            dY2 = 0;
        else
            dY2 = Yscale * V(Structure.DOF(2, Nodes(ElmNd(j,2)).ID), IX(i));
        end
        Ydef(2,j) = Nodes(ElmNd(j,2)).Ycoord + dY2;
    end
    line (Xdef,Ydef,'Color',[1,0,0],'LineWidth',[0.5]); % plot deformed elements
    NdXdisp = zeros(1, Structure.NumNodes); NdYdisp = zeros(1, Structure.NumNodes);
    for k = 1:Structure.NumNodes
        if Structure.DOF(1, Nodes(k).ID) == -1
            dX = 0;
        else
            dX = Xscale * V(Structure.DOF(1, Nodes(k).ID), IX(i));
        end
        NdXdisp(1,k) = Nodes(k).Xcoord + dX;
        if Structure.DOF(2, Nodes(k).ID) == -1
            dY = 0;
        else
            dY = Yscale * V(Structure.DOF(2, Nodes(k).ID), IX(i));
        end
        NdYdisp(1,k) = Nodes(k).Ycoord + dY;
    end
    scatter (NdXdisp,NdYdisp,MarkerSize,'s','MarkerEdgeColor','r'); % plot displaced nodes
    axis equal
    box on
    title(['Mode No.',num2str(i),' [Period = ',num2str(Structure.System_periods(i),'%6.4f' ),' sec]']);
end
