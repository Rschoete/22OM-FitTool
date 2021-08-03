%!!!!!!!!!!!!!!!!!!!Move file to main directory before run!!!!!!!!!!!!!!!!!!!!!!
close all
clear all
clc

nr= '200409';
solverMethod = 'ms';
saveorrun = 'save';  % save: Save settings below to input or start: start fit procedure

Input.filename = sprintf('Target4SBsuro_%s.mat',nr);
Input.TargetValues_sp = [15:49];
load(fullfile('./Targets',Input.filename))
Input.TargetValues_tr = find([Target.TauRecov(:).Vm]<0&[Target.TauRecov(:).Il]>=100);
settings = {};
addpath(genpath('./Functions'))


%changes to flags

%changes to standard Input fit
figpos = [60,187,1791,757];
powera = 1;
powerb = 1;
settings = horzcat(settings,{'powera',powera,'powerb',powerb,'figpos',figpos});
% changes standard optimization otpions
maxIter = 10e5;
maxTime = 24*3600;
DisplayIter = 'iter';
OdeOpts = {'RelTol',1e-8,'AbsTol',1e-8,'Maxstep',100e-6};
settings = horzcat(settings,{'maxIter',maxIter,'DisplayIter',DisplayIter,'OdeOpts',OdeOpts,'maxTime',maxTime});

%changes Tau
options_TauO.msnr = 200;
options_TauDA.msnr = 200;
settings = horzcat(settings,{'options_TauO',options_TauO,'Options_TauDA',options_TauDA});
%changes fit GODA
wGODA.type = 'equal';
options_GODA.tend = 1.5;
options_GODA.dt = 1.5e-4;

options_GODA.msnr = 300;
options_GODA.solverMethod = solverMethod;

options_GODA.psoBC.options.Generations = 100;
%options_GODA.psoBC.options.PopulationSize = 700;       %default is 10*nvars
options_GODA.psoBC.options.PlotFcns = {@psoplotbestf,@psoplotswarmsurf} ;
options_GODA.psoBC.options.TimeLimit = maxTime;

options_GODA.ga.options.Generations = 100;
%options_GODA.ga.options.PopulationSize = 300;      %default is 10*nvars
options_GODA.ga.options.PlotFcns = {@gaplotbestf,@gaplotstopping};
options_GODA.ga.options.TimeLimit = maxTime;

settings = horzcat(settings,{'options_GODA',options_GODA,'wGODA',wGODA});
% changes final optimization
fopt_solverMethod = solverMethod;

options_fopt.ftype = 'frac_initb';
options_fopt.fLB = 1/10;
options_fopt.fUB = 1/10;
options_fopt.FB.wtype = 'input';
options_fopt.FB.w = repmat([1000,1000,1000,1,1,10],length(Input.TargetValues_sp),1);
options_fopt.FB.wtr = repmat(10,length(Input.TargetValues_tr),1);
options_fopt.ms.msnr = 700;
options_fopt.FB.nrr_featextr = 4;
options_fopt.FB.tend_sp = 2;
options_fopt.Outfun.plot_flag = 0;

options_fopt.psoBC.options.Generations = 100;
%options_fopt.psoBC.options.PopulationSize = 700;       %default is 10*nvars
options_fopt.psoBC.options.PlotFcns = {@psoplotbestf,@psoplotswarmsurf} ;
options_fopt.psoBC.options.TimeLimit = maxTime;

options_fopt.ga.options.Generations = 100;
%options_fopt.ga.options.PopulationSize = 700;      %default is 10*nvars
options_fopt.ga.options.PlotFcns = {@gaplotbestf,@gaplotstopping};
options_fopt.ga.options.TimeLimit = maxTime;

settings = horzcat(settings,{'options_fopt',options_fopt,'fopt_solverMethod',fopt_solverMethod});


fun.Gr_Williams.fun=@(p,V) (10.6408-14.6408*exp(-(V-p(1))/42.7671));
fun.Gr_Williams.nrpara=1;
fun.Gr_Williams.LB=[-100];
fun.Gr_Williams.X0=[0];
fun.Gr_Williams.UB=[100];
fun.Gr_Williams.var=1;

OinfIfun_names = {'symLogisticslog'};
TauOIfun_names = {'symLogisticsloginv'};
DAinfIfun_names = {'symLogisticslogmina'};
TauDAIfun_names = {'doublesymLogisticslogmina'};
TauOVfun_names = {'Logistics'};
TauDAVfun_names = {'Logistics'};
Grectfun_names = {'Gr_Williams'};
combTauO = @(funV,funI) funV.*funI;
combTauDA = @(funV,funI) funV.*funI;

settings = horzcat(settings,{'newfun',fun,'OinfI',OinfIfun_names,'TauOI',TauOIfun_names,'DAinfI',DAinfIfun_names,...
    'TauDAI',TauDAIfun_names,'TauOV',TauOVfun_names,'TauDAV',TauDAVfun_names,'Grect',Grectfun_names,'combTauO',combTauO,'combTauDA',combTauDA});


Savename = fullfile('./Inputs/H134Rsuro/',sprintf('H134Rsuro_s%s_t%s_d%s_nw.mat',solverMethod,nr,datestr(now,'yymmddHH')));

switch saveorrun
case 'save'
    save(Savename,'Input','settings')
case 'run'
    out = fit22HH(Input,settings);
end
