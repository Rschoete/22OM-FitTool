%create input file used by fit tool
close all
clear all
clc

test_flag = 1; % when active reduces time limits and number of particles and or starting points
saveorrun = 'save';  % save: Save settings below to input or start: start fit procedure

% Example target H134R
nr= '200414';
Input.filename = sprintf('TargetH134R_%s.mat',nr);
load(fullfile('./Targets',Input.filename));
settings = {};
addpath(genpath('./Functions'))

% Example of some hyperparameters which can be adjusted
% good practice to add a string which indicates changes made to fit
% procedure. Do not forget to add str to suffix_ex at line 150

solverMethods = {'psoBC'}; %'ms','ga',
%wGODAtypes = {'equal','input','input'};
% wGODAw = {[1,1,1],[1,2,5],[1,1,2]};
% wGODAstr = {'','_iwl','_iws'};
wGODAtypes = {'input'};
wGODAw = {[1,2,5]};
wGODAstr = {'_iwl'};
DAlimflag = [0,1];
DAlimstr = {'','_DAlim'};
if test_flag
    Timelim = [10*60];%,12*3600];
    Timelim_str = {'_Time10m'};%,'_Time12h'};
else
    Timelim = [12*3600];%,12*3600];
    Timelim_str = {'_Time12h'};%,'_Time12h'};
end

CombTauO_str = {'_prod','_recsum'};
%CombTauDA_str = {'_prod'},'_recsum'};
CombTauO = {@(funV,funI) funV.*funI, @(funV,funI) ((funV).^(-1)+(funI).^(-1)).^(-1)};
CombTauDA = {@(funV,funI) funV.*funI,@(funV,funI) ((funV).^(-1)+(funI).^(-1)).^(-1)};
%CombTauDA = {@(funV,funI) ((funV).^(-1)+(funI).^(-1)).^(-1)};
CombTauDA_str = {'_prod','_recsum'};


% loop over different hyperparmeters.
% for more info on parameters see declare_parameters function under
% Functions/Subfunction
for is = 1:length(solverMethods)
    for iwg = 1:length(wGODAtypes)
        for idl = 1:length(DAlimflag)
            for iTl = 1:length(Timelim)
                for iCTO = 1:length(CombTauO)
                    for iCTDA = 1
                    iCTDA = iCTO;

solverMethod = solverMethods{is};





%changes to flags

%changes to standard Input fit
figpos = [60,187,1791,757];
powera = 1;
powerb = 1;
settings = horzcat(settings,{'powera',powera,'powerb',powerb,'figpos',figpos});

% changes standard optimization otpions
if test_flag
    maxIter = 1e2;
else
    maxIter = 1e6;
end
maxTime = Timelim(iTl);%24*3600;
DisplayIter = 'iter';
StartPointsToRun = 'all';
OdeOpts = {'RelTol',1e-8,'AbsTol',1e-8,'Maxstep',100e-6};
settings = horzcat(settings,{'maxIter',maxIter,'DisplayIter',DisplayIter,'OdeOpts',OdeOpts,'maxTime',maxTime,'StartPointsToRun',StartPointsToRun});

%changes Tau
if test_flag
    options_TauO.msnr = 20;
    options_TauDA.msnr = 20;
else
    options_TauO.msnr = 2000;
    options_TauDA.msnr = 2000;
end
settings = horzcat(settings,{'options_TauO',options_TauO,'Options_TauDA',options_TauDA});
%changes fit GODA
inclIratio = 1;
wGODA.type = wGODAtypes{iwg};
wGODA.w = repmat(wGODAw{iwg},length(Target.SingleP),1);
options_GODA = struct();
options_GODA.tend = 1.5;
options_GODA.dt = 1.5e-4;


options_GODA.solverMethod ='ms';
if test_flag
    options_GODA.msnr = 100;
else
    options_GODA.msnr = 3000;
end
options_GODA.extrac_nonlcon.IO_bigger_flag = 1;
options_GODA.extrac_nonlcon.Olim_bigger = 0.60;
options_GODA.extrac_nonlcon.IO_bigger = 5500;
if DAlimflag(idl)
    options_GODA.extrac_nonlcon.IDA_bigger_flag = 1;
options_GODA.extrac_nonlcon.DAlim_bigger = 0.6;
options_GODA.extrac_nonlcon.IDA_bigger = 10;
end

if test_flag
    options_GODA.psoBC.options.Generations = 10;
    options_GODA.psoBC.options.PopulationSize = 30;       %default is 10*nvars
else
    options_GODA.psoBC.options.Generations = 1000;
    options_GODA.psoBC.options.PopulationSize = 300;       %default is 10*nvars
end
%options_GODA.psoBC.options.PlotFcns = {@psoplotbestf,@psoplotswarmsurf} ;
options_GODA.psoBC.options.TimeLimit = maxTime;
if test_flag
    options_GODA.ga.options.Generations = 10;
    options_GODA.ga.options.PopulationSize = 30;      %default is 10*nvars
else
    options_GODA.ga.options.Generations = 1000;
    options_GODA.ga.options.PopulationSize = 300;      %default is 10*nvars
end
%options_GODA.ga.options.PlotFcns = {@gaplotbestf,@gaplotstopping};
options_GODA.ga.options.TimeLimit = maxTime;

settings = horzcat(settings,{'options_GODA',options_GODA,'wGODA',wGODA,'inclIratio',inclIratio});



% changes final optimization
fopt_solverMethod = solverMethod;

options_fopt.ftype = 'frac_initb';
options_fopt.fLB = 1/10;
options_fopt.fUB = 1/10;
options_fopt.FB.wtype = 'input';
options_fopt.FB.w = repmat([1000,1000,1000,10*wGODAw{iwg}],length(Target.SingleP),1);%repmat([1000,1000,1000,1,1,10],length(Target.SingleP),1);%
options_fopt.FB.wtr = repmat(20,length(Target.TauRecov),1);
if test_flag
    options_fopt.ms.msnr = 10;
else
    options_fopt.ms.msnr = 1000;
end
options_fopt.FB.nrr_featextr = 4;
options_fopt.FB.tend_sp = 2;
options_fopt.Outfun.plot_flag = 0;
options_fopt.extrac_nonlcon.IO_bigger_flag = 1;
options_fopt.extrac_nonlcon.Olim_bigger = 0.60;
options_fopt.extrac_nonlcon.IO_bigger = 5500;
options_fopt.FB.method_tr = 'currenttraces';

if test_flag
    options_fopt.psoBC.options.Generations = 10;
else
    options_fopt.psoBC.options.Generations = 1000;
end
%options_fopt.psoBC.options.PopulationSize = 10;       %default is 10*nvars
%options_fopt.psoBC.options.PlotFcns = {@psoplotbestf}%,@psoplotswarmsurf} ;
options_fopt.psoBC.options.TimeLimit = maxTime;

if test_flag
    options_fopt.ga.options.Generations = 10;
else
    options_fopt.ga.options.Generations = 1000;
end
%options_fopt.ga.options.PopulationSize = 700;      %default is 10*nvars
%options_fopt.ga.options.PlotFcns = {@gaplotbestf,@gaplotstopping};
options_fopt.ga.options.TimeLimit = maxTime;

settings = horzcat(settings,{'options_fopt',options_fopt,'fopt_solverMethod',fopt_solverMethod});


OinfIfun_names = {'symLogisticslog'};
TauOIfun_names = {'symLogisticsloginv'};
DAinfIfun_names = {'symLogisticslogmina'};
TauDAIfun_names = {'doublesymLogisticslogmina7'};
TauOVfun_names = {'Logistics'};
TauDAVfun_names = {'Logistics'};
Grectfun_names = {'Grect'};
combTauO = CombTauO{iCTO};%@(funV,funI) funV.*funI;
combTauDA = CombTauDA{iCTDA};%@(funV,funI) funV.*funI;

settings = horzcat(settings,{'OinfI',OinfIfun_names,'TauOI',TauOIfun_names,'DAinfI',DAinfIfun_names,...
    'TauDAI',TauDAIfun_names,'TauOV',TauOVfun_names,'TauDAV',TauDAVfun_names,'Grect',Grectfun_names,'combTauO',combTauO,'combTauDA',combTauDA});


suffix_ex = [DAlimstr{idl},wGODAstr{iwg},Timelim_str{iTl},CombTauO_str{iCTO},CombTauDA_str{iCTDA},'_TDA_dllma7'];
if test_flag
    suffix_ex = [suffix_ex,'_test'];
end

% (out)comment below what you want to do. Currently only input file
% created fit is not started

Savename = fullfile('./Inputs/H134Rexp/',sprintf('H134R_s%s_t%s_d%s%s.mat',solverMethod,nr,datestr(now,'yymmddHH'),suffix_ex));
switch saveorrun
    case 'save'
        save(Savename,'Input','settings')
    case 'run'
        out = fit22HH(Input,settings);
end
                    end
                end
            end
        end
    end
end