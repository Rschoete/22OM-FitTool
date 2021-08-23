%Compile script
%create standalone application of fit tool. Usefull when only one matlab
%license available. This way it can be run without license after standalone application is created.
%License is necessary prior compilation

%paths below schould be adjusted according to own system and folders
addpath('$user\PSO_BC');% necessary Particle Swarm optimization with boundary conditions


% -I includes folder path and compiles matlab functions within, does not
% include other files => use -a. -a increases size of standalone
% application -I is sufficient and input files can still be changed if put
% into Inputs, no recompilation is necessary while this would be the case
% with -a
mcc -mv runfit22HH_LIS.m -I ./Functions -I ./Functions/Subfunctions -I ./Inputs -I ./Targets -d ./Compilations