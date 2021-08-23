%!!!!!!!!!!!!!!!!!!!Move file to main directory before run!!!!!!!!!!!!!!!!!!!!!!
close all
clear all
clc

saveorrun = 'save';  % save: Save settings below to input or start: start fit procedure

solverMethods = {'ms'}; %'ms','ga',
%wGODAtypes = {'equal','input','input'};
% wGODAw = {[1,1,1],[1,2,5],[1,1,2]};
% wGODAstr = {'','_iwl','_iws'};
wGODAtypes = {'equal'};
wGODAw = {[1,1,1]};
wGODAstr = {''};
DAlimflag = [0];
DAlimstr = {''};

for is = 1:length(solverMethods)
    for iwg = 1:length(wGODAtypes)
        for idl = 1:length(DAlimflag)

solverMethod = solverMethods{is};
nr= '200414';

Input.filename = sprintf('TargetH134R_%s.mat',nr);
load(fullfile('./Targets',Input.filename));
settings = {};
addpath(genpath('./Functions'))

Vvals = unique([Target.SingleP(:).Vm]);

%changes to flags

%changes to standard Input fit
figpos = [60,187,1791,757];
powera = 1;
powerb = 1;
settings = horzcat(settings,{'powera',powera,'powerb',powerb,'figpos',figpos});
% changes standard optimization otpions
maxIter = 1000;
maxTime = 24*3600;
DisplayIter = 'iter';
StartPointsToRun = 'all';
OdeOpts = {'RelTol',1e-8,'AbsTol',1e-8,'Maxstep',100e-6};
settings = horzcat(settings,{'maxIter',maxIter,'DisplayIter',DisplayIter,'OdeOpts',OdeOpts,'maxTime',maxTime,'StartPointsToRun',StartPointsToRun});

%changes Tau
options_TauO.msnr = 2000;
options_TauDA.msnr = 2000;
settings = horzcat(settings,{'options_TauO',options_TauO,'Options_TauDA',options_TauDA});
%changes fit GODA
inclIratio = 1;
wGODA.type = wGODAtypes{iwg};
wGODA.w = repmat(wGODAw{iwg},length(Target.SingleP),1);
options_GODA = struct();
options_GODA.tend = 1.5;
options_GODA.dt = 1.5e-4;


options_GODA.solverMethod ='ms';

options_GODA.msnr = 3000;
options_GODA.extrac_nonlcon.IO_bigger_flag = 1;
options_GODA.extrac_nonlcon.Olim_bigger = 0.60;
options_GODA.extrac_nonlcon.IO_bigger = [Vvals(:),ones(size(Vvals(:)))*5500];
if DAlimflag(idl)
    options_GODA.extrac_nonlcon.IDA_bigger_flag = 1;
options_GODA.extrac_nonlcon.DAlim_bigger = 0.6;
options_GODA.extrac_nonlcon.IDA_bigger = [Vvals(:),ones(size(Vvals(:)))*10];
end

options_GODA.psoBC.options.Generations = 1000;
options_GODA.psoBC.options.PopulationSize = 300;       %default is 10*nvars
%options_GODA.psoBC.options.PlotFcns = {@psoplotbestf,@psoplotswarmsurf} ;
options_GODA.psoBC.options.TimeLimit = maxTime;

options_GODA.ga.options.Generations = 1000;
options_GODA.ga.options.PopulationSize = 300;      %default is 10*nvars
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
options_fopt.ms.msnr = 10;
options_fopt.FB.nrr_featextr = 4;
options_fopt.FB.tend_sp = 2;
options_fopt.Outfun.plot_flag = 0;
options_fopt.extrac_nonlcon.IO_bigger_flag = 1;
options_fopt.extrac_nonlcon.Olim_bigger = 0.60;
options_fopt.extrac_nonlcon.IO_bigger = [Vvals(:),ones(size(Vvals(:)))*5500];
options_fopt.FB.method_tr = 'currenttraces';

options_fopt.psoBC.options.Generations = 1000;
%options_fopt.psoBC.options.PopulationSize = 10;       %default is 10*nvars
%options_fopt.psoBC.options.PlotFcns = {@psoplotbestf}%,@psoplotswarmsurf} ;
options_fopt.psoBC.options.TimeLimit = maxTime;

options_fopt.ga.options.Generations = 100;
%options_fopt.ga.options.PopulationSize = 700;      %default is 10*nvars
%options_fopt.ga.options.PlotFcns = {@gaplotbestf,@gaplotstopping};
options_fopt.ga.options.TimeLimit = maxTime;

settings = horzcat(settings,{'options_fopt',options_fopt,'fopt_solverMethod',fopt_solverMethod});


% OinfIfun_names = {'symLogisticslog'};
% TauOIfun_names = {'symLogisticsloginv'};
% DAinfIfun_names = {'symLogisticslogmina'};
% TauDAIfun_names = {'doublesymLogisticslogmina'};
% TauOVfun_names = {'Logistics'};
% TauDAVfun_names = {'Logistics'};
OinfIfun_names = {'OinfIV'};
TauOIfun_names = {'TauOIV'};
DAinfIfun_names = {'DAinfdIlIV'};
TauDAIfun_names = {'TauDAdIlIV'};
TauOVfun_names = {'One'};
TauDAVfun_names = {'One'};
Grectfun_names = {'Grect'};
combTauO = @(funV,funI) funV.*funI;
combTauDA = @(funV,funI) funV.*funI;

settings = horzcat(settings,{'OinfI',OinfIfun_names,'TauOI',TauOIfun_names,'DAinfI',DAinfIfun_names,...
    'TauDAI',TauDAIfun_names,'TauOV',TauOVfun_names,'TauDAV',TauDAVfun_names,'Grect',Grectfun_names,'combTauO',combTauO,'combTauDA',combTauDA});


suffix_ex = [DAlimstr{idl},wGODAstr{iwg}];

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

%%



%changes to standard Input fit
powera = 1;
powerb = 1;
settings = horzcat(settings,{'powera',powera,'powerb',powerb});
% changes standard optimization otpions
maxIter = 1e5;
maxTime = 12*3600;
DisplayIter = 'iter';
OdeOpts = {'RelTol',1e-8,'AbsTol',1e-8,'Maxstep',100e-6};
settings = horzcat(settings,{'maxIter',maxIter,'DisplayIter',DisplayIter,'OdeOpts',OdeOpts,'maxTime',maxTime});

%changes Tau
options_TauO.msnr = 200;
options_TauDA.msnr = 200;
settings = horzcat(settings,{'options_TauO',options_TauO,'Options_TauDA',options_TauDA});
%changes fit GODA
options_GODA.tend = 1.5;
options_GODA.dt = 1.5e-4;
options_GODA.msnr = 300;
settings = horzcat(settings,{'options_GODA',options_GODA});
% changes final optimization
options_fopt.FB.wtype = 'rel';
options_fopt.ms.msnr = 700;
options_fopt.FB.nrr_featextr = 4;

settings = horzcat(settings,{'options_fopt',options_fopt});


OinfIfun_names = {'symLogisticslog'};
TauOIfun_names = {'symLogisticsloginv'};
DAinfIfun_names = {'symLogisticslogmina'};
TauDAIfun_names = {'doublesymLogisticslogmina'};
TauOVfun_names = {'Logistics'};
TauDAVfun_names = {'Logistics'};
Grectfun_names = {'Grect'};
combTauO = @(funV,funI) funV.*funI;
combTauDA = @(funV,funI) ((funV).^(-1)+(funI).^(-1)).^(-1);

settings = horzcat(settings,{'OinfI',OinfIfun_names,'TauOI',TauOIfun_names,'DAinfI',DAinfIfun_names,...
    'TauDAI',TauDAIfun_names,'TauOV',TauOVfun_names,'TauDAV',TauDAVfun_names,'Grect',Grectfun_names,'combTauO',combTauO,'combTauDA',combTauDA});

%out(7) = fit22HH(Input,settings);
