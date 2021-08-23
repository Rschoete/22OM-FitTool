%run22fitMM
%!!!!!!!!!!!!!!!!!!!Move file to main directory before run!!!!!!!!!!!!!!!!!!!!!!
clear all
close all
clc

saveorrun = 'save';  % save: Save settings below to input or start: start fit procedure

solverMethods = {'psoBC','ms'};
wton = [10];
wton_str = {'_tOnL'};
fitOs = {'tO-tDA-GODA','GODA-tDA-tO'};
fitOs_str = {'_tOtDAGODA','_GODAtDAtO'};
CombTauO_str = {'_prod','_recsum'};
%CombTauDA_str = {'_prod'},'_recsum'};
CombTauO = {@(funV,funI) funV.*funI, @(funV,funI) ((funV).^(-1)+(funI).^(-1)).^(-1)};
%CombTauDA = {@(funV,funI) funV.*funI},@(funV,funI) ((funV).^(-1)+(funI).^(-1)).^(-1)};
CombTauDA = {@(funV,funI) ((funV).^(-1)+(funI).^(-1)).^(-1)};
CombTauDA_str = {'_recsum'};

for is = 1%:length(solverMethods)
    for iwton = 1%:length(wton)
        for ifOs = 1%:length(fitOs)
            for iCTO = 1:length(CombTauO)
                for iCTDA = 1:length(CombTauDA)

nr= '200416_modToff2';
solverMethod = solverMethods{is};

Input.filename = sprintf('TargetMM%s.mat',nr);
Input.TargetValues_sp = [1:7,15:19];
settings = {};
addpath(genpath('./Functions'))
%changes to flags

%changes to standard Input fit
figpos = [60,187,1791,757];
powera = 1;
powerb = 1;
settings = horzcat(settings,{'powera',powera,'powerb',powerb,'figpos',figpos});
% changes standard optimization otpions
maxIter = 1e6;
maxTime = 8*3600;
DisplayIter = 'iter';
OdeOpts = {'RelTol',1e-8,'AbsTol',1e-8,'Maxstep',100e-6};
settings = horzcat(settings,{'maxIter',maxIter,'DisplayIter',DisplayIter,'OdeOpts',OdeOpts,'maxTime',maxTime});
% change fit order

fitOrder = fitOs{ifOs};
settings = horzcat(settings,{'fitOrder',fitOrder});

%changes Tau
options_TauO.msnr = 2000;
options_TauDA.msnr = 2000;
settings = horzcat(settings,{'options_TauO',options_TauO,'Options_TauDA',options_TauDA});

%changes fit GODA
load(fullfile('./Targets',Input.filename))
wGODA.type = 'input';
wGODA.w = repmat([1,50],length(Input.TargetValues_sp),1);
options_GODA = struct();
options_GODA.tend = 1.5;
options_GODA.dt = 1.5e-4;

options_GODA.msnr = 3000;
options_GODA.solverMethod = 'ms';

options_GODA.psoBC.options.Generations = 1000;
options_GODA.psoBC.options.PopulationSize = 300;       %default is 10*nvars
%options_GODA.psoBC.options.PlotFcns = {@psoplotbestf,@psoplotswarmsurf} ;
options_GODA.psoBC.options.TimeLimit = maxTime;

options_GODA.ga.options.Generations = 1000;
options_GODA.ga.options.PopulationSize = 300;      %default is 10*nvars
%options_GODA.ga.options.PlotFcns = {@gaplotbestf,@gaplotstopping};
options_GODA.ga.options.TimeLimit = maxTime;

options_GODA.extrac_nonlcon.IO_smaller_flag = 1;
options_GODA.extrac_nonlcon.IO_smaller = [4000];
options_GODA.extrac_nonlcon.Olim_smaller = [0.5];
options_GODA.extrac_nonlcon.IO_bigger_flag = 0;
options_GODA.extrac_nonlcon.IO_bigger = [];
options_GODA.extrac_nonlcon.Olim_bigger = [];
settings = horzcat(settings,{'options_GODA',options_GODA,'wGODA',wGODA});
% changes final optimization
%options_fopt.FB.wtype = 'sd';
fopt_solverMethod = solverMethod;

options_fopt.ftype = 'frac_initb';
options_fopt.fLB = 1/10;
options_fopt.fUB = 1/10;
options_fopt.FB.wtype = 'input';
options_fopt.FB.w = repmat([1000*wton(iwton),1000,1000,1/25,1,250],length(Input.TargetValues_sp),1);
options_fopt.FB.wtr = repmat(10,length(Target.TauRecov),1);
options_fopt.ms.msnr = 1000;
options_fopt.FB.nrr_featextr = 4;
options_fopt.extrac_nonlcon.IO_smaller_flag = 1;
options_fopt.extrac_nonlcon.IO_smaller = [4000];
options_fopt.extrac_nonlcon.Olim_smaller = [0.5];
options_fopt.extrac_nonlcon.IO_bigger_flag = 0;
options_fopt.extrac_nonlcon.IO_bigger = [];
options_fopt.extrac_nonlcon.Olim_bigger = [];
options_fopt.Outfun.plot_flag = 0;


options_fopt.psoBC.options.Generations = 100;
%options_fopt.psoBC.options.PopulationSize = 700;       %default is 10*nvars
%options_fopt.psoBC.options.PlotFcns = {@psoplotbestf,@psoplotswarmsurf} ;
options_fopt.psoBC.options.TimeLimit = maxTime;

options_fopt.ga.options.Generations = 100;
%options_fopt.ga.options.PopulationSize = 700;      %default is 10*nvars
%options_fopt.ga.options.PlotFcns = {@gaplotbestf,@gaplotstopping};
options_fopt.ga.options.TimeLimit = maxTime;

settings = horzcat(settings,{'options_fopt',options_fopt,'fopt_solverMethod',fopt_solverMethod});

% extra fun
fun.Logistics_TI.fun=@(p,V) p(1)./(1+exp(-(V-p(2))./p(3)))+p(4);
fun.Logistics_TI.nrpara=4;
fun.Logistics_TI.LB=[0 -100 -100,0];
fun.Logistics_TI.X0=[4e-3 0 20 7e-3];
fun.Logistics_TI.UB=[5e-2,100,100,5e-2];
fun.Logistics_TI.var=1;

fun.Logistics_TO.fun=@(p,V) p(1)./(1+exp(-(V-p(2))./p(3)))+p(4);
fun.Logistics_TO.nrpara=4;
fun.Logistics_TO.LB=[0 -100 -100,0];
fun.Logistics_TO.X0=[0.5e-3 0 30 0.9e-3];
fun.Logistics_TO.UB=[1e-2,100,100,1e-2];
fun.Logistics_TO.var=1;

fun.symLogisticsloginv_TO.fun=@(p,Il) (p(3)./(p(1).*Il.^(1./(p(2)*log(10)))+1));
fun.symLogisticsloginv_TO.nrpara=3;
fun.symLogisticsloginv_TO.LB=[0,0,0];
fun.symLogisticsloginv_TO.X0=[0.6,0.8,0.05];
fun.symLogisticsloginv_TO.UB=[2,2,0.5];
fun.symLogisticsloginv_TO.var=2;

fun.symLogisticsloginv_TI.fun=@(p,Il) (p(3)./(p(1).*Il.^(1./(p(2)*log(10)))+1));
fun.symLogisticsloginv_TI.nrpara=3;
fun.symLogisticsloginv_TI.LB=[0,0,0];
fun.symLogisticsloginv_TI.X0=[0.6,0.5,7*1.2708];
fun.symLogisticsloginv_TI.UB=[2,2,20];
fun.symLogisticsloginv_TI.var=2;

fun.symLogisticslog_O.fun=@(p,I) (1./(1+p(1).*I.^(-1./(p(2)*log(10)))));
fun.symLogisticslog_O.LB=[100 0];
fun.symLogisticslog_O.X0=[10000,0.4];
fun.symLogisticslog_O.UB=[3e4,10];
fun.symLogisticslog_O.nrpara=2;
fun.symLogisticslog_O.var=2;

fun.symLogisticsloga_DA.fun=@(p,I) 1-(p(3)./(1+p(1).*I.^(-1./(p(2)*log(10)))));
fun.symLogisticsloga_DA.LB=[10,0,0];
fun.symLogisticsloga_DA.X0=[10000,0.4,0.05];
fun.symLogisticsloga_DA.UB=[3e4,10,0.5];
fun.symLogisticsloga_DA.nrpara=3;
fun.symLogisticsloga_DA.var=2;

fun.Gr_cte.fun=@(p,V) p(1)*(V-p(2));
fun.Gr_cte.nrpara=2;
fun.Gr_cte.LB=[0,-100];
fun.Gr_cte.X0=[30,0];
fun.Gr_cte.UB=[300,100];
fun.Gr_cte.var=1;

fun.symLogisticslogmina.fun=@(p,I) 1-(p(3)./(1+exp(p(1)/p(2)).*I.^(-1./(p(2)*log(10)))));
fun.symLogisticslogmina.LB=[-100,0,0.80];
fun.symLogisticslogmina.X0=[1,1,0.9];
fun.symLogisticslogmina.UB=[10,20,1];
fun.symLogisticslogmina.nrpara=3;
fun.symLogisticslogmina.var=2;


% OinfIfun_names = {'symLogisticslog','symLogisticslog_O'};
% TauOIfun_names = {'symLogisticsloginv','symLogisticsloginv_TO'};
% DAinfIfun_names = {'symLogisticslogmina','symLogisticsloga_DA'};
% TauDAIfun_names = {'doublesymLogisticslogmina','symLogisticsloginv_TI'};
% TauOVfun_names = {'Logistics','Logistics_TO'};
% TauDAVfun_names = {'Logistics','Logistics_TI'};
% Grectfun_names = {'GHH','Gr_cte'};

OinfIfun_names = {'symLogisticslog'};
TauOIfun_names = {'symLogisticsloginv'};
DAinfIfun_names = {'symLogisticslogmina'};
TauDAIfun_names = {'doublesymLogisticslogmina'};
TauOVfun_names = {'Logistics'};
TauDAVfun_names = {'Logistics'};
Grectfun_names = {'GHH'};

% OinfIfun_names = {'symLogisticslog_O'};
% TauOIfun_names = {'symLogisticsloginv_TO'};
% DAinfIfun_names = {'symLogisticsloga_DA'};
% TauDAIfun_names = {'symLogisticsloginv_TI'};
% TauOVfun_names = {'Logistics_TO'};
% TauDAVfun_names = {'Logistics_TI'};
% Grectfun_names = {'Gr_cte'};


combTauO = CombTauO{iCTO};%@(funV,funI) funV.*funI;
combTauDA = CombTauDA{iCTDA};%@(funV,funI) funV.*funI;

settings = horzcat(settings,{'newfun',fun,'OinfI',OinfIfun_names,'TauOI',TauOIfun_names,'DAinfI',DAinfIfun_names,...
    'TauDAI',TauDAIfun_names,'TauOV',TauOVfun_names,'TauDAV',TauDAVfun_names,'Grect',Grectfun_names,'combTauO',combTauO,'combTauDA',combTauDA});

suffix_ex = [wton_str{iwton},fitOs_str{ifOs},CombTauO_str{iCTO},CombTauDA_str{iCTDA}];

Savename = fullfile('./Inputs/MM/',sprintf('MM_s%s_t%s_d%s%s.mat',solverMethod ,nr,datestr(now,'yymmddHH'),suffix_ex));


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
