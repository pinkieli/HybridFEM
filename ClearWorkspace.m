clc;close all;clrok=input('[FEM] Clear Workspace First? y/[n]: ', 's');
if ~isempty(clrok)&&((clrok=='y'||clrok=='Y'))
    clear;
end;
global TARGET RUNMODE;