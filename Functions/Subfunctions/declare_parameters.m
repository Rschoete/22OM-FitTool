function [checkDouble_flag,display_flag,errorbar_flag,memorySave_flag,plotFeat_flag,plotTau_flag,plot_flag,plotTauRecov_flag,...
    TauRecov_flag,sort_flag,time_flag,inclIratio_flag,corrTauDAIl0_flag,powera,powerb,OdeOpts,featNames,featNames_SD,featNamesTr,featNamesTr_SD,...
    featNamesT_plot,featNamesI_plot,IVdepend_targets,recal_opt_FBIVdepend,fitOrder,options_TauO,options_TauDA,...
    wGODA,options_GODA,fopt_select,fopt_objMethod,fopt_solverMethod,options_fopt,Colors,type,fig_pos] =...
    declare_parameters(Target_sp,Target_tr,varargin)
tf_string = {'false','true'};

% default settings
useParallel = 1; %'off', 'always'
StartPointsToRun = 'all'; %'bounds','bounds-ineqs'
% flags
checkDouble_flag = 1;
display_flag = 1;
errorbar_flag = 1;
HPC_flag = 0;
TauRecov_flag = 1;
plotTauRecov_flag = 1;
plotFeat_flag = 1;
plotTau_flag = 1;
plot_flag = 1;
memorySave_flag = 0;
sort_flag = 1;
time_flag = 1;
inclIratio_flag = 0;
corrTauDAIl0_flag = 1; % correction of TauDA for Il = 0.
%TauDAIl0 estimate of Tau recov however only valid if DAinf ~=0 otherwise
%correction needed Taurec ~= TauDA0(1-ln(1/(1-Iratio)) ==> if TauDA0 set equal
%to Taurec real Tau rec is underestimated.
%real relationship exp(-Trec/TDA0) =
%(1-A*B-exp(-1))/(A*(DAon-1)*exp(-dtp2/TDA))*(DAinf+(1-DAinf)*exp(-dtp1/TDA))
%relationship above is obtained if under assumption dtp2 = dtp1 DAon =
%DAinf and dtp1>TO see paper https://doi.org/10.1101/2020.11.10.376939


% Do not change following flags!!
redef_msopt_flag = 0;
redef_optimopt_flag = 0;
redef_OdeOpts_flag = 0;
adjFields_roptFBIV_flag = 0;
adjFields_optTauO_flag = 0;
adjFields_optTauDA_flag = 0;
adjFields_optGODA_flag = 0;
adjFields_optfopt_flag = 0;

% Standard Input fit
powera = 1;
powerb = 1;

% standard optimization options
maxIter = 1000;
maxTime = inf;
maxFunEval = 1e6;
DisplayIter = 'final';             %'final','iter','off'
funTol = 1e-10;
XTol = 1e-10;
OdeOpts = {'RelTol',1e-12,'AbsTol',1e-12,'Maxstep',25e-6};

% used in CheckDouble
featNames = {'Ipeak','Iss','Imean','Iratio','Imean','TauOn','TauOff','TauInact','ttp','I'};
featNames_SD = cellfun(@(str) [str,'_SD'], featNames, 'UniformOutput', false);
featNamesTr = {'TauRecov'};
featNamesTr_SD = {'TauRecov_SD'};

% feature names to plot in featplot
featNamesT_plot = {'TauOff','TauOn','TauInact'};
featNamesI_plot = {'Ipeak','Iss','Iratio'};

% Fit TauO and TauDA

IVdepend_targets = 'features';     %feautres, currenttraces, postcurrenttraces

% feature based: recal TAuOn Off Inact rec to correct powera and b
recal_opt_FBIVdepend.TauO.optimoptions = {'lsqcurvefit','Algorithm','trust-region-reflective',...
    'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
recal_opt_FBIVdepend.TauO.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};
recal_opt_FBIVdepend.TauO.msnr = 20;
% optional options
%recal_opt_FBIVdepend.TauO.ub = 1;  recal_opt_FBIVdepend.TauO.lb = 0;
%recal_opt_FBIVdepend.TauO.tend = 0.5;  recal_opt_FBIVdepend.TauO.dt = 1e-6;

recal_opt_FBIVdepend.TauDA.optimoptions = {'lsqcurvefit','Algorithm','trust-region-reflective',...
    'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter'};
recal_opt_FBIVdepend.TauDA.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};
recal_opt_FBIVdepend.TauDA.msnr = 20;
% optional options
%recal_opt_FBIVdepend.TauDA.lb = 0; recal_opt_FBIVdepend.TauDA.lbr = 0;
%recal_opt_FBIVdepend.TauDA.ub = 1; recal_opt_FBIVdepend.TauDA.ubr = 30;
%recal_opt_FBIVdepend.TauDA.tend = 0.5; recal_opt_FBIVdepend.TauDA.tendr = 15;
%recal_opt_FBIVdepend.TauDA.dt = 1e-2; recal_opt_FBIVdepend.TauDA.dtr = 1e-4;


% fit order
fitOrder = 'tO-tDA-GODA';           %tO-tDA-GODA,  GODA-tDA-tO, noFirstFit

% options fitTauO
options_TauO.msnr = 200;                %100 per function TaO comb two
options_TauO.solverMethod = 'lsqcurvefit'; % lsqcurvefit','fmincon'
options_TauO.optimoptions = {'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
options_TauO.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};

% options fitTauDA
options_TauDA.msnr = 200;               % 100 per function
options_TauDA.solverMethod = 'lsqcurvefit'; % lsqcurvefit','fmincon'
options_TauDA.optimoptions = {'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
options_TauDA.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};

% options fitGODA
wGODA.type = 'sd';          %'sd', 'input', 'equal','rel'
wGODA.w = 1;
options_GODA.I_nonlcon = []; %e.g., logspace(0,4,5)
options_GODA.V1_nonlcon = []; % e.g., [10,20,40]
options_GODA.V2_nonlcon = []; %e.g,[-80,-60,-40,-20,-10]
options_GODA.extrac_nonlcon.IO_smaller_flag = 0;
options_GODA.extrac_nonlcon.IO_smaller = [];
options_GODA.extrac_nonlcon.Olim_smaller = [];
options_GODA.extrac_nonlcon.IO_bigger_flag = 0;
options_GODA.extrac_nonlcon.IO_bigger = [];
options_GODA.extrac_nonlcon.Olim_bigger = [];
options_GODA.extrac_nonlcon.IDA_smaller_flag = 0;
options_GODA.extrac_nonlcon.IDA_smaller = [];
options_GODA.extrac_nonlcon.DAlim_smaller = [];
options_GODA.extrac_nonlcon.IDA_bigger_flag = 0;
options_GODA.extrac_nonlcon.IDA_bigger = [];
options_GODA.extrac_nonlcon.DAlim_bigger = [];
options_GODA.extraceq_nonlcon.VG = [];
options_GODA.tend = 2;
options_GODA.dt = 2e-4;

options_GODA.solverMethod = 'ms';       %'ms' 'psobc' 'ga'
% GODA fit solverMethod;
options_GODA.msnr = 300;
options_GODA.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};
options_GODA.optimoptions = {'fmincon','Algorithm','sqp',...
    'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
% GODA: PSO with nonlinBounds
options_GODA.psoBC.Aineq = []; options_GODA.psoBC.Aeq = [];
options_GODA.psoBC.bineq = [];  options_GODA.psoBC.beq = [];
options_GODA.psoBC.options.UseParallel = useParallel;
options_GODA.psoBC.options.Generations = 1000;
options_GODA.psoBC.options.Display = DisplayIter;
options_GODA.psoBC.options.TolFun = funTol;
%options_GODA.psoBC.options.OutputFcns = @stopPSO;

% GODA: GA
options_GODA.ga.Aineq = []; options_GODA.ga.Aeq = [];
options_GODA.ga.bineq = [];  options_GODA.ga.beq = [];
options_GODA.ga.options.UseParallel = useParallel;
options_GODA.ga.options.Generations = 1000;
options_GODA.ga.options.Display = DisplayIter;
options_GODA.ga.options.TolFun = funTol;
%options_GODA.ga.options.OutputFcns = @stopPSO;  %should be applicable here as well

options_GODA.hyboptions.UseParallel = useParallel;
options_GODA.hyboptions.Display = DisplayIter;
options_GODA.hyboptions.MaxIter = maxIter;
options_GODA.hyboptions.FunctionTolerance = funTol;
options_GODA.hyboptions.OptimalityTolerance = funTol;
options_GODA.hyboptions.StepTolerance = XTol;
options_GODA.hyboptions.XTolerance = XTol;
options_GODA.hyboptions.MaxFunctionEvaluations = maxFunEval;
options_GODA.hyboptions.MaxIterations = maxIter;
options_GODA.hyboptions.Algorithm = 'sqp';
options_GODA.hybfun = 'fmincon';


% options final optimization
fopt_select = 'all';               %'all' or 'best'
fopt_objMethod = 'features';        %'features' or 'currenttraces'
fopt_solverMethod = 'ms';           %'ms', 'psobc','pso' 'ga';
options_fopt.ftype = 'scale';       %'scale', 'add', 'frac_initb'
options_fopt.fLB = 10^(-1/2);
options_fopt.fUB = 10^(1/2);
options_fopt.incl_tr_flag = 1;
options_fopt.CT.method_CT = 'normalized'; %'normalized','zscore','weighted','normalized+weighted'
options_fopt.CT.obj = 'rms';                %rms or q2
% Input required!
%   options_fopt.CT.intervals
%otptional
%   options_fopt.CT.w
%   options_fopt.CT.wtr

options_fopt.FB.method_tr = 'currenttraces';        %estimpowerb','noestimpowerb',currenttraces' (first to tend to underestimate if Iratio >0
options_fopt.FB.trestim.optimoptions = {'lsqcurvefit','Algorithm','trust-region-reflective',...
    'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
options_fopt.FB.trestim.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};
options_fopt.FB.trestim.msnr = 20;
options_fopt.FB.wtype = 'sd';   %'sd','rel','input'
options_fopt.FB.dt_tr = 1e-3;  %1ms
options_fopt.FB.dxx    = 1e-3; %1ms
options_fopt.FB.dt = 1e-4; %1ms
options_fopt.FB.nrr_featextr = 20;
options_fopt.FB.plot_featextr = 0;
options_fopt.I_nonlcon = []; %[logspace(0,4,5)]
options_fopt.V1_nonlcon = []; %[10,20,40]
options_fopt.V2_nonlcon = []; %[-80,-60,-40,-20,-10]
options_fopt.extrac_nonlcon.IO_smaller_flag = 0;
options_fopt.extrac_nonlcon.IO_smaller = [];
options_fopt.extrac_nonlcon.Olim_smaller = [];
options_fopt.extrac_nonlcon.IO_bigger_flag = 0;
options_fopt.extrac_nonlcon.IO_bigger = [];
options_fopt.extrac_nonlcon.Olim_bigger = [];
options_fopt.extrac_nonlcon.IDA_smaller_flag = 0;
options_fopt.extrac_nonlcon.IDA_smaller = [];
options_fopt.extrac_nonlcon.DAlim_smaller = [];
options_fopt.extrac_nonlcon.IDA_bigger_flag = 0;
options_fopt.extrac_nonlcon.IDA_bigger = [];
options_fopt.extrac_nonlcon.DAlim_bigger = [];
options_fopt.extraceq_nonlcon.VG = [];
options_fopt.Outfun.OdeOpts = OdeOpts;
options_fopt.Outfun.plot_flag = 0;

% fopt_solverMethod: MultiStart
options_fopt.ms.optimoptions = {'fmincon','Algorithm','sqp',...
    'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
options_fopt.ms.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};
options_fopt.ms.msnr = 700;

% fopt_solverMethod: PSO with nonlinBounds
options_fopt.psoBC.Aineq = []; options_fopt.psoBC.Aeq = [];
options_fopt.psoBC.bineq = [];  options_fopt.psoBC.beq = [];
options_fopt.psoBC.options.UseParallel = useParallel;
options_fopt.psoBC.options.Generations = 1000;
options_fopt.psoBC.options.Display = DisplayIter;
options_fopt.psoBC.options.TolFun = funTol;
%options_fopt.psoBC.options.OutputFcns = @stopPSO;

% fopt_solverMethod: GA
options_fopt.ga.Aineq = []; options_fopt.ga.Aeq = [];
options_fopt.ga.bineq = [];  options_fopt.ga.beq = [];
options_fopt.ga.options.UseParallel = useParallel;
options_fopt.ga.options.Generations = 1000;
options_fopt.ga.options.Display = DisplayIter;
options_fopt.ga.options.TolFun = funTol;
%options_fopt.ga.options.OutputFcns = @stopPSO;  %should be applicable here as well

options_fopt.hyboptions.UseParallel = useParallel;
options_fopt.hyboptions.Display = DisplayIter;
options_fopt.hyboptions.MaxIter = maxIter;
options_fopt.hyboptions.FunctionTolerance = funTol;
options_fopt.hyboptions.OptimalityTolerance = funTol;
options_fopt.hyboptions.StepTolerance = XTol;
options_fopt.hyboptions.XTolerance = XTol;
options_fopt.hyboptions.MaxFunctionEvaluations = maxFunEval;
options_fopt.hyboptions.MaxIterations = maxIter;
options_fopt.hyboptions.Algorithm = 'sqp';
options_fopt.hybfun = 'fmincon';




%Input required
%    options_fopt.FB.intervals
%    options_fopt.FB.w TauOn,TauInact,TAuOff,Ipeak, Iss,Iratio
%   options_fopt.FB.wtr
%optional
%   options_fopt.FB.tend_sp
%   options_fopt.FB.tend_tr
%   options_fopt.FB.trestim.lb
%   options_fopt.FB.trestim.ub
%   options_fopt.FB.trestim.dt
%   options_fopt.FB.trestim.tend

% figures
Colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840; 1 ,0 ,0; 0, 1, 0; 0, 0, 1];
type = {'-','o-','--','-.','-*'};
%fig_pos = [-1603,206,1404,727];


% Change standard input parameters
if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'HPC'))
            HPC_flag = varargin{find(strcmpi(varargin,'HPC'))+1};
            if HPC_flag
                display_flag = 0;
                plot_flag = 0;
                plotFeat_flag = 0;
                plotTau_flag = 0;
                plotTauRecov_flag = 0;
                DisplayIter = 'final';
                maxIter = 1e6;
                redef_msopt_flag = 1;
                redef_optimopt_flag = 1;
            end
        end
        if any(strcmpi(varargin,'plot'))
            plot_flag = varargin{find(strcmi(varargin,'plot'))+1};
            if plot_flag
                plotFeat_flag = 1;
                plotTau_flag = 1;
                plotTauRecov_flag = 1;
            end
        end
        if any(strcmpi(varargin,'TauRecov'))
            TauRecov_flag = varargin{find(strcmpi(varargin,'TauRecov'))+1};
            if ~TauRecov_flag
                plotTauRecov_flag = 0;
                options_fopt.incl_tr_flag = 0;
            end
        end
        if any(strcmpi(varargin,'checkDouble'))
            checkDouble_flag = varargin{find(strcmpi(varargin,'checkDouble'))+1};
        end
        if any(strcmpi(varargin,'Colors'))
            Colors = varargin{find(strcmpi(varargin,'Colors'))+1};
        end
        if any(strcmpi(varargin,'corrTauDAIl0_flag'))
            corrTauDAIl0_flag = varargin{find(strcmpi(varargin,'corrTauDAIl0_flag'))+1};
        end
        if any(strcmpi(varargin,'Display'))
            display_flag = varargin{find(strcmpi(varargin,'Display'))+1};
        end
        if any(strcmpi(varargin,'DisplayIter'))
            DisplayIter = varargin{find(strcmpi(varargin,'DisplayIter'))+1};
            redef_msopt_flag = 1;
        end
        if any(strcmpi(varargin,'errorbar'))
            errorbar_flag = varargin{find(strcmpi(varargin,'errorbar'))+1};
        end
        if any(strcmpi(varargin,'featNames'))
            featNames = varargin{find(strcmpi(varargin,'featNames'))+1};
            featNames_SD = cellfun(@(str) [str,'_SD'], featNames, 'UniformOutput', false);
        end
        if any(strcmpi(varargin,'featNames_tr'))
            featNamesTr = varargin{find(strcmpi(varargin,'featNames_tr'))+1};
            featNamesTr_SD = cellfun(@(str) [str,'_SD'], featNamesTr, 'UniformOutput', false);
        end
        if any(strcmpi(varargin,'featNamesI'))
            featNamesI_plot = varargin{find(strcmpi(varargin,'featNamesI'))+1};
        end
        if any(strcmpi(varargin,'featNamesT'))
            featNamesT_plot = varargin{find(strcmpi(varargin,'featNamesT'))+1};
        end
        if any(strcmpi(varargin,'figpos'))
            fig_pos = varargin{find(strcmpi(varargin,'figpos'))+1};
        end
        if any(strcmpi(varargin,'fitOrder'))
            fitOrder = varargin{find(strcmpi(varargin,'fitOrder'))+1};
        end
        if any(strcmpi(varargin,'fopt_objMethod'))
            fopt_objMethod = varargin{find(strcmpi(varargin,'fopt_objMethod'))+1};
        end
        if any(strcmpi(varargin,'fopt_select'))
            fopt_select = varargin{find(strcmpi(varargin,'fopt_select'))+1};
        end
        if any(strcmpi(varargin,'fopt_solverMethod'))
            fopt_solverMethod = varargin{find(strcmpi(varargin,'fopt_solverMethod'))+1};
        end
        if any(strcmpi(varargin,'funTol'))
            funTol = varargin{find(strcmpi(varargin,'funTol'))+1};
            redef_msopt_flag = 1;
        end
        if any(strcmpi(varargin,'IVdepend_targets'))
            IVdepend_targets = varargin{find(strcmpi(varargin,'IVdepend_targets'))+1};
        end
        if any(strcmpi(varargin,'inclIratio'))
            inclIratio_flag = varargin{find(strcmpi(varargin,'inclIratio'))+1};
        end
        if any(strcmpi(varargin,'maxIter'))
            maxIter = varargin{find(strcmpi(varargin,'maxIter'))+1};
            redef_optimopt_flag = 1;
        end
        if any(strcmpi(varargin,'maxFunEval'))
            maxFunEval = varargin{find(strcmpi(varargin,'maxFunEval'))+1};
            redef_optimopt_flag = 1;
        end
        if any(strcmpi(varargin,'maxTime'))
            maxTime = varargin{find(strcmpi(varargin,'maxTime'))+1};
            redef_msopt_flag = 1;
        end
        if any(strcmpi(varargin,'memorySave'))
            memorySave_flag = varargin{find(strcmpi(varargin,'memorySave'))+1};
        end
        if any(strcmpi(varargin,'OdeOpts'))
            OdeOpts = varargin{find(strcmpi(varargin,'OdeOpts'))+1};
            redef_OdeOpts_flag = 1;
        end
        if any(strcmpi(varargin,'options_TauO'))
            options_TauO_input = varargin{find(strcmpi(varargin,'options_TauO'))+1};
            adjFields_optTauO_flag = 1;
        end
        if any(strcmpi(varargin,'options_TauDA'))
            options_TauDA_input = varargin{find(strcmpi(varargin,'options_TauDA'))+1};
            adjFields_optTauDA_flag = 1;
        end
        if any(strcmpi(varargin,'options_GODA'))
            options_GODA_input = varargin{find(strcmpi(varargin,'options_GODA'))+1};
            adjFields_optGODA_flag = 1;
        end
        if any(strcmpi(varargin,'options_fopt'))
            options_fopt_input = varargin{find(strcmpi(varargin,'options_fopt'))+1};
            adjFields_optfopt_flag = 1;
        end
        if any(strcmpi(varargin,'plotFeat'))
            plotFeat_flag = varargin{find(strcmpi(varargin,'plotFeat'))+1};
        end
        if any(strcmpi(varargin,'plotTau'))
            plotTau_flag = varargin{find(strcmi(varargin,'plotTau'))+1};
        end
        if any(strcmpi(varargin,'plotTauRecov'))
            plotTauRecov_flag = varargin{find(strcmpi(varargin,'plotTauRecov'))+1};
        end
        if any(strcmpi(varargin,'powera'))
            powera = varargin{find(strcmpi(varargin,'powera'))+1};
        end
        if any(strcmpi(varargin,'powerb'))
            powerb = varargin{find(strcmpi(varargin,'powerb'))+1};
        end
        if any(strcmpi(varargin,'recal_opt'))
            recal_opt_input = varargin{find(strcmpi(varargin,'recal_opt_FBIVdepend'))+1};
            adjFields_roptFBIV_flag = 1;
        end
        if any(strcmpi(varargin,'Sort'))
            sort_flag = varargin{find(strcmpi(varargin,'Sort'))+1};
        end
        if any(strcmpi(varargin,'StartPointsToRun'))
            StartPointsToRun = varargin{find(strcmpi(varargin,'StartPointsToRun'))+1};
            redef_msopt_flag = 1;
        end
        if any(strcmpi(varargin,'time_flag'))
            time_flag = varargin{find(strcmpi(varargin,'time_flag'))+1};
        end
        if any(strcmpi(varargin,'linetype'))
            type = varargin{find(strcmpi(varargin,'linetype'))+1};
        end
        if any(strcmpi(varargin,'wGODA'))
            wGODA_new = varargin{find(strcmpi(varargin,'wGODA'))+1};
            wGODA = adjFields(wGODA,wGODA_new,1);
        end
        if any(strcmpi(varargin,'XTol'))
            XTol = varargin{find(strcmpi(varargin,'XTol'))+1};
            redef_msopt_flag = 1;
        end
    end
end



if redef_msopt_flag
    recal_opt_FBIVdepend.TauO.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};
    recal_opt_FBIVdepend.TauDA.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};
    options_TauO.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};
    options_TauDA.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};
    options_GODA.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};
    options_fopt.FB.trestim.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};
    options_fopt.ms.msopt = {'Display',DisplayIter,'UseParallel',useParallel,'FunctionTolerance',funTol,'XTolerance',XTol,'MaxTime',maxTime,'StartPointsToRun',StartPointsToRun};

end
if redef_optimopt_flag
    recal_opt_FBIVdepend.TauO.optimoptions = {'lsqcurvefit','Algorithm','trust-region-reflective',...
        'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
    recal_opt_FBIVdepend.TauDA.optimoptions = {'lsqcurvefit','Algorithm','trust-region-reflective',...
        'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter'};
    options_TauO.optimoptions = {'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
    options_TauDA.optimoptions = {'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
    options_GODA.optimoptions = {'fmincon','Algorithm','sqp',...
        'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
    options_fopt.ms.optimoptions = {'fmincon','Algorithm','sqp',...
        'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter,'Display','off'};
    options_fopt.FB.trestim.optimoptions = {'lsqcurvefit','Algorithm','trust-region-reflective',...
        'MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
end
if redef_msopt_flag || redef_optimopt_flag
    options_GODA.hyboptions.UseParallel = useParallel;
    options_GODA.hyboptions.Display = DisplayIter;
    options_GODA.hyboptions.MaxIter = maxIter;
    options_GODA.hyboptions.FunctionTolerance = funTol;
    options_GODA.hyboptions.OptimalityTolerance = funTol;
    options_GODA.hyboptions.StepTolerance = XTol;
    options_GODA.hyboptions.XTolerance = XTol;
    options_GODA.hyboptions.MaxFunctionEvaluations = maxFunEval;
    options_GODA.hyboptions.MaxIterations = maxIter;

    options_GODA.psoBC.options.Display = DisplayIter;
    options_GODA.psoBC.options.TolFun = funTol;

    options_GODA.ga.options.Display = DisplayIter;
    options_GODA.ga.options.TolFun = funTol;

    options_fopt.hyboptions.UseParallel = useParallel;
    options_fopt.hyboptions.Display = DisplayIter;
    options_fopt.hyboptions.MaxIter = maxIter;
    options_fopt.hyboptions.FunctionTolerance = funTol;
    options_fopt.hyboptions.OptimalityTolerance = funTol;
    options_fopt.hyboptions.StepTolerance = XTol;
    options_fopt.hyboptions.XTolerance = XTol;
    options_fopt.hyboptions.MaxFunctionEvaluations = maxFunEval;
    options_fopt.hyboptions.MaxIterations = maxIter;

    options_fopt.psoBC.options.Display = DisplayIter;
    options_fopt.psoBC.options.TolFun = funTol;

    options_fopt.ga.options.Display = DisplayIter;
    options_fopt.ga.options.TolFun = funTol;
end
if redef_OdeOpts_flag
    options_fopt.Outfun.OdeOpts = OdeOpts;
end
if adjFields_roptFBIV_flag
    recal_opt_FBIVdepend = adjFields(recal_opt_FBIVdepend,recal_opt_input,2);
end
if adjFields_optTauO_flag
    options_TauO = adjFields(options_TauO,options_TauO_input,1);
end
if adjFields_optTauDA_flag
    options_TauDA = adjFields(options_TauDA,options_TauDA_input,1);
end
if adjFields_optGODA_flag
    options_GODA = adjFields(options_GODA,options_GODA_input,2);
end
if adjFields_optfopt_flag
    options_fopt = adjFields(options_fopt,options_fopt_input,3);
end

if strcmpi(options_GODA.solverMethod,'psoBC')
    hybridopt = optimoptions(options_GODA.hybfun);
    Hfns=fieldnames(hybridopt);
    Hfns_n = fieldnames(options_GODA.hyboptions);
    for i=1:length(Hfns)
        if any(strcmpi(Hfns{i},Hfns_n))
            fieldname = Hfns_n{strcmpi(Hfns{i},Hfns_n)};
            hybridopt.(Hfns{i})=options_GODA.hyboptions.(fieldname);
        end
    end
    options_GODA.psoBC.options.HybridFcn = {str2func(options_GODA.hybfun),hybridopt};
    options_GODA = rmfield(options_GODA,{'hybfun','hyboptions'});
end

if strcmpi(options_GODA.solverMethod,'ga')
    hybridopt = optimoptions(options_GODA.hybfun);
    Hfns=fieldnames(hybridopt);
    Hfns_n = fieldnames(options_GODA.hyboptions);
    for i=1:length(Hfns)
        if any(strcmpi(Hfns{i},Hfns_n))
            fieldname = Hfns_n{strcmpi(Hfns{i},Hfns_n)};
            hybridopt.(Hfns{i})=options_GODA.hyboptions.(fieldname);
        end
    end
    options_GODA.ga.options.HybridFcn = {str2func(options_GODA.hybfun),hybridopt};
    options_GODA = rmfield(options_GODA,{'hybfun','hyboptions'});
end

if strcmpi(fopt_solverMethod,'psoBC')
    hybridopt = optimoptions(options_fopt.hybfun);
    Hfns=fieldnames(hybridopt);
    Hfns_n = fieldnames(options_fopt.hyboptions);
    for i=1:length(Hfns)
        if any(strcmpi(Hfns{i},Hfns_n))
            fieldname = Hfns_n{strcmpi(Hfns{i},Hfns_n)};
            hybridopt.(Hfns{i})=options_fopt.hyboptions.(fieldname);
        end
    end
    options_fopt.psoBC.options.HybridFcn = {str2func(options_fopt.hybfun),hybridopt};
    options_fopt = rmfield(options_fopt,{'hybfun','hyboptions'});
end

if strcmpi(fopt_solverMethod,'ga')
    hybridopt = optimoptions(options_fopt.hybfun);
    Hfns=fieldnames(hybridopt);
    Hfns_n = fieldnames(options_fopt.hyboptions);
    for i=1:length(Hfns)
        if any(strcmpi(Hfns{i},Hfns_n))
            fieldname = Hfns_n{strcmpi(Hfns{i},Hfns_n)};
            hybridopt.(Hfns{i})=options_fopt.hyboptions.(fieldname);
        end
    end
    options_fopt.ga.options.HybridFcn = {str2func(options_fopt.hybfun),hybridopt};
    options_fopt = rmfield(options_fopt,{'hybfun','hyboptions'});
end


%adjust errorbar based on Target
%Check if TauRecov data included
if TauRecov_flag
    fnms_Targettr = fieldnames(Target_tr);
    targetTrSD_flag = ~any(contains(fnms_Targettr,'_SD','IgnoreCase',true));
else
    %check if Target_tr is empty ==> should be otherwise trow warning
    if ~isempty(Target_tr)
        warning('Target_tr is not empty but not TauRecov_Flag = 0')
    end
    targetTrSD_flag = 1;
end
fnms_Targetsp = fieldnames(Target_sp);
if ~any(contains(fnms_Targetsp,'_SD','IgnoreCase',true)) || targetTrSD_flag
    if errorbar_flag
        warning('errorbar_flag selected but no SD in target')
        pause(0.1)
        fprintf('\nerrorbar_flag off now\n')
        errorbar_flag = 0;
    end
end



if display_flag
    fprintf('\nflags\n------\n')
    fprintf('chekcDoubleflag: %s\n',tf_string{checkDouble_flag+1})
    fprintf('errorbarflag: %s \n',tf_string{errorbar_flag+1});
    fprintf('HPC flag; %s\n',tf_string{HPC_flag+1});
    fprintf('plotflag: %s \n',tf_string{plot_flag+1});
    fprintf('plot features: %s\n',tf_string{plotFeat_flag+1});
    fprintf('plot Tau: %s\n',tf_string{plotTau_flag+1});
    fprintf('sortflag: %s \n',tf_string{sort_flag+1});
    fprintf('memorySave flag: %s\n',tf_string{memorySave_flag+1});
    fprintf('timeflag: %s \n',tf_string{time_flag+1});
    fprintf('include TauRecovery data: %s\n',tf_string{TauRecov_flag+1})
    fprintf('plot Tau Recov as well: %s \n',tf_string{plotTauRecov_flag+1});
    fprintf('Outfun plot flag: %s \n',tf_string{options_fopt.Outfun.plot_flag+1})
    fprintf('inlcude tau Recov in final optimization: %s \n',tf_string{options_fopt.incl_tr_flag+1});
    fprintf('include Iratio in goda: %s \n',tf_string{inclIratio_flag+1});
    fprintf('correct TauDA0 not being equal to Taurec: %s \n',tf_string{corrTauDAIl0_flag+1});
    fprintf('\n adj_flags\n---------\n')
    fprintf('adj Fields recal options FB IV dependence: %s\n',tf_string{adjFields_roptFBIV_flag+1})
    fprintf('adj Fields options TauO: %s\n',tf_string{adjFields_optTauO_flag+1})
    fprintf('adj Fields options TauDA: %s\n',tf_string{adjFields_optTauDA_flag+1})
    fprintf('adj Fields options G Oinf DAin: %s\n',tf_string{adjFields_optGODA_flag+1})
    fprintf('adj Fields options final optimisation: %s\n',tf_string{adjFields_optfopt_flag+1})
    fprintf('\n fit options\n----------\n')
    fprintf('maximum iterations: %5.2e\n',maxIter)
    fprintf('maximum time: %5.2e\n',maxTime)
    fprintf('maximum function evaluations: %5.2e\n',maxFunEval)
    fprintf('Start points ms: %s\n',StartPointsToRun);
    fprintf('function tolerance: %5.2e\n',funTol)
    fprintf('X tolerance: %5.2e\n',XTol)
    fprintf('nr multistart recal options FBIVdepend TauO: %i\n',recal_opt_FBIVdepend.TauO.msnr)
    fprintf('nr multistart recal options FBIVdepend TauDA: %i\n',recal_opt_FBIVdepend.TauDA.msnr)
    fprintf('nr multistart options TauO: %i\n',options_TauDA.msnr)
    fprintf('nr multistart options TauDA: %i\n',options_TauDA.msnr)
    fprintf('nr multistart options G Oinf DAin: %i\n',options_GODA.msnr)
    fprintf('nr multistart final optimization FB ms: %i\n',options_fopt.ms.msnr)
    fprintf('nr multistart final optimization FB TauRecov estimation: %i\n',options_fopt.FB.trestim.msnr)
    fprintf('\nMethods and inputs\n----------\n')
    fprintf('power a: %i \n',powera);
    fprintf('power b: %i \n',powerb);
    fprintf('featurenames included: %s\n',strjoin(featNames))
    fprintf('Taurecovery featnames: %s\n',strjoin(featNamesTr))
    fprintf('Target type IVdependencies: %s\n',IVdepend_targets)
    fprintf('fit Order: %s\n',fitOrder);
    fprintf('solver method TauO: %s\n',options_TauO.solverMethod)
    fprintf('solver method TauDA: %s\n',options_TauDA.solverMethod)
    fprintf('Options LB UB scaling: %s\n',options_fopt.ftype)
    fprintf('Ode Options: %s = %5.2e\n',OdeOpts{:})
    fprintf('weights G Oinf DAinf: %s\n',wGODA.type)
    fprintf('Solver Method GODA: %s\n',options_GODA.solverMethod);
    fprintf('final optimization selection: %s\n',fopt_select)
    fprintf('final optimization objective method: %s\n',fopt_objMethod)
    fprintf('final optimization solverMethod: %s\n',fopt_solverMethod)
    if strcmpi(fopt_objMethod,'currenttraces')
        fprintf('method currenttraces: %s\n',options_fopt.CT.method_CT)
        fprintf('costfunciton: %s\n',options_fopt.CT.obj)

    elseif strcmpi(fopt_objMethod,'features')
        fprintf('TauRecov method: %s\n',options_fopt.FB.method_tr)
        fprintf('weights fopt FB: %s\n',options_fopt.FB.wtype)
    end
end

% final checkups
if strcmpi(wGODA.type,'input')
    if inclIratio_flag
        if size(wGODA.w,1)~=length(Target_sp) || size(wGODA.w,2)~=3
            error('weights wGoda incorrect!')
        end
    else
        if size(wGODA.w,1)~=length(Target_sp) || size(wGODA.w,2)~=2
            error('weights wGoda incorrect!')
        end
    end
end
if strcmpi(fopt_objMethod,'currenttraces')
    if options_fopt.incl_tr_flag
        if ~isfield(options_fopt.CT,'intervals')
            try
                options_fopt.CT.intervals = Target_tr(1).Intervals;
            catch
                error('no interval field')
            end
        end
    end
    if contains(options_fopt.CT.method_CT,'weighted')
        if isfield(options_fopt.CT,'w')
            if size(options_fopt_CT.w,1)~=length(Target_sp)
                error('size weights doenst match')
            end
        else
            error('no weights included')
        end
        if isfield(options_fopt.CT,'wtr')
            if size(options_fopt_CT.wtr,1)~=length(Target_tr)
                error('size weights doesnt match')
            end
        else
            error('no weights tr included')
        end
    end
elseif strcmpi(fopt_objMethod,'features')
    switch options_fopt.FB.wtype
        case 'sd'
            options_fopt.FB.w = 1./[[Target_sp(:).TauOn_SD]',[Target_sp(:).TauInact_SD]',[Target_sp(:).TauOff_SD]',...
                [Target_sp(:).Ipeak_SD]',[Target_sp(:).Iss_SD]',[Target_sp(:).Iratio_SD]'];
        case 'rel'
            options_fopt.FB.w = 1./[[Target_sp(:).TauOn]',[Target_sp(:).TauInact]',[Target_sp(:).TauOff]',...
                [Target_sp(:).Ipeak]',[Target_sp(:).Iss]',[Target_sp(:).Iratio]'];
        case 'input'
            %included above under adjFields_optfopt_flag
            if ~adjFields_optfopt_flag
                error('no input weights')
            end
        case 'rel&input'
            options_fopt.FB.w = options_fopt.FB./[[Target_sp(:).TauOn]',[Target_sp(:).TauInact]',[Target_sp(:).TauOff]',...
                [Target_sp(:).Ipeak]',[Target_sp(:).Iss]',[Target_sp(:).Iratio]'];
        case 'sd&input'
            options_fopt.FB.w = options_fopt.FB.w./[[Target_sp(:).TauOn_SD]',[Target_sp(:).TauInact_SD]',[Target_sp(:).TauOff_SD]',...
                [Target_sp(:).Ipeak_SD]',[Target_sp(:).Iss_SD]',[Target_sp(:).Iratio_SD]'];
        otherwise
            error('false w type')
    end
    if isfield(options_fopt.FB,'w')
        if size(options_fopt.FB.w,1)~=length(Target_sp) && size(options_fopt.FB.w,2)~=6
            error('size weights doenst match')
        end
    else
        error('no weights included')
    end
    if options_fopt.incl_tr_flag
        switch options_fopt.FB.wtype
            case 'sd'
                options_fopt.FB.wtr = [Target_tr(:).TauRecov_SD]';
            case 'rel'
                options_fopt.FB.wtr = [Target_tr(:).TauRecov]';
            case 'input'
                %included above under adjFields_optfopt_flag
                if ~adjFields_optfopt_flag
                    error('no input weights')
                end
            otherwise
                error('false w type')
        end
        if isfield(options_fopt.FB,'wtr')
            if size(options_fopt.FB.wtr,1)~=length(Target_tr)
                error('size weights doesnt match')
            end
        else
            error('no weights tr included')
        end

        if ~isfield(options_fopt.FB,'intervals')
            try
                options_fopt.FB.intervals = Target_tr(1).Intervals;
            catch
                try
                options_fopt.FB.intervals = Target_tr(1).intervals;
                catch
                error('no interval field')
                end
            end
        end
    end
else
    error('false fopt_objMethods')
end

Colors = @(idx) Colors(mod(idx-1,size(Colors,1))+1,:);
type = @(idx) type{mod(idx-1,length(type))+1};

end