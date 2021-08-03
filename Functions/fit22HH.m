function Out = fit22HH(Input,varargin)
% This functton fits a HH type model to Target data of an opsin.
% Input should be a structure that contains the filename under
% Input.filename. The target file to load should be stored under the
% Targets folder in the current working directory. The Target.mat file
% should contain two fields, Singlep and TauRecov, which contain the data
% of single pulse and two pulse experiments respectively. The input
% structure can contain also another field which allow for selection of a
% subgroup of data to fit to. These should be stored in TargetValues_sp and
% TargetValues_tr for single pulse and TauRecov respectively.
%
% varargin contains all alternative settings and parameters according to
% the stringname,indicative for the parameter to change, and value principle.
% Varargin is thus always even

starttimer = tic;
%unpack varargin (possible to give a single cell to varargin containing the info as described above) 
if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end

% select target values
try
    load(fullfile('./Targets/',Input.filename))
catch
    try
        load(Input.filename)

    catch
        error('invalid target file')
    end
end
Out.Targetfilename = Input.filename;
Out.Input = Input;

Target_sp = Target.SingleP;
Target_tr = Target.TauRecov;
if isfield(Input,'TargetValues_sp')
    Target_sp = Target_sp(Input.TargetValues_sp);
    Out.TargetValues_sp = Input.TargetValues_sp;
else
    fprintf('\n all target values are included in single pulse \n')
end
if isfield(Input,'TargetValues_tr')
    Target_tr = Target_tr(Input.TargetValues_tr);
    Out.TargetValues_tr = Input.TargetValues_tr;
else
    fprintf('\n all target values are included in tau rec \n')
end

[checkDouble_flag,display_flag,errorbar_flag,memorySave_flag,plotFeat_flag,plotTau_flag,plot_flag,plotTauRecov_flag,...
    TauRecov_flag,sort_flag,time_flag,inclIratio_flag,corrTauDAIl0_flag,powera,powerb,OdeOpts,featNames,featNames_SD,featNamesTr,featNamesTr_SD,...
    featNamesT_plot,featNamesI_plot,IVdepend_targets,recal_opt_FBIVdepend,fitOrder,options_TauO,options_TauDA,...
    wGODA,options_GODA,fopt_select,fopt_objMethod,fopt_solverMethod,options_fopt,Colors,type,fig_pos] =...
    declare_parameters(Target_sp,Target_tr,varargin);


if time_flag
    Timer(1).value = toc(starttimer);
    Timer(1).comment = 'loaded parameters';
end

%Different Il and Vm0 values
Il_valsp = unique([Target_sp(:).Il]);
Vm_valsp = unique([Target_sp(:).Vm]);
if TauRecov_flag
    Il_valtr = unique([Target_tr(:).Il]);
    Vm_valtr = unique([Target_tr(:).Vm]);
else
    Il_valtr = [];
    Vm_valtr = [];
end

%Check for doubles in Target file
if checkDouble_flag
    % single pulse
    Targetnew = Target_sp;
    indextocut = [];
    fnTn = fieldnames(Targetnew);
    for iVm = Vm_valsp
        sindex = find([Target_sp(:).Vm]==iVm);
        sTarget = Target_sp(sindex);
        [n, bin] = histc([sTarget(:).Il], unique([sTarget(:).Il]));
        multiple = find(n==2);
        if any(n>2)
            error('more than two equal samples adjust code')
        end
        if ~isempty(multiple)
            index = find(ismember(bin, multiple));
            sTarget = Target_sp(sindex(index));
            indextocut = horzcat(indextocut,sindex(index(2)));
            for ifN = 1:length(featNames)
                if any(strcmpi(fnTn,featNames{ifN}))
                    if strcmpi(featNames{ifN},'I') && length(sTarget(1).(featNames{ifN}))~=length(sTarget(2).(featNames{ifN}))
                        idx = (sTarget(2).nsamples>sTarget(1).nsamples)+1;
                        Targetnew(sindex(index(1))).(featNames{ifN}) = sTarget(idx).(featNames{ifN});
                        Targetnew(sindex(index(1))).nsamplesI = sTarget(idx).nsamples;
                        Targetnew(sindex(index(1))).time = sTarget(idx).time;
                    else
                    Targetnew(sindex(index(1))).(featNames{ifN}) = 1/sum([sTarget(:).nsamples])*...
                        (sTarget(1).nsamples*sTarget(1).(featNames{ifN}) + sTarget(2).nsamples*sTarget(2).(featNames{ifN}));
                    if any(strcmpi(fnTn,featNames_SD{ifN}))
                        Targetnew(sindex(index(1))).(featNames_SD{ifN}) = 1/sum([sTarget(:).nsamples])*...
                            (sTarget(1).nsamples*(sTarget(1).(featNames_SD{ifN})^2+sTarget(1).(featNames{ifN})^2)+...
                            sTarget(2).nsamples*(sTarget(2).(featNames_SD{ifN})^2+sTarget(2).(featNames{ifN})^2))-...
                            Targetnew(sindex(index(1))).(featNames{ifN})^2;
                        Targetnew(sindex(index(1))).(featNames_SD{ifN}) = sqrt(Targetnew(sindex(index(1))).(featNames_SD{ifN}));
                    else
                        warning('no SD info')
                    end
                    end
                else
                    warning('not all fieldnames present')
                end
            end
            Targetnew(sindex(index(1))).('nsamples') = sTarget(1).nsamples+sTarget(2).nsamples;
        end
    end
    Targetnew(indextocut) = [];
    Target_sp = Targetnew;
    clearvars('Targetnew','sTarget','index','multiple','n','bin','sindex')
    %adjust weights as well
    if strcmpi(wGODA.type,'input')
        wGODA.w(indextocut,:) = [];
    end
    switch fopt_objMethod
        case 'currenttraces'
            if contains(options_fopt.CT.method_CT,'weighted','IgnoreCase',true)
                options_fopt.CT.w(indextocut,:) = [];
            end
        case 'features'
            options_fopt.FB.w(indextocut,:) = [];
        otherwise; error('why?');
    end
    if TauRecov_flag
        %TauRecov
        Targetnew = Target_tr;
        indextocut = [];
        fnTn = fieldnames(Targetnew);
        for iVm = Vm_valtr
            sindex = find([Target_tr(:).Vm]==iVm);
            sTarget = Target_tr(sindex);
            [n, bin] = histc([sTarget(:).Il], unique([sTarget(:).Il]));
            multiple = find(n==2);
            if any(n>2)
                error('more than two equal samples adjust code')
            end
            if ~isempty(multiple)
                index = find(ismember(bin, multiple));
                sTarget = Target_tr(sindex(index));
                indextocut = horzcat(indextocut,sindex(index(2)));
                for ifN = 1:length(featNamesTr)
                    if any(strcmpi(fnTn,featNamesTr{ifN}))
                        Targetnew(sindex(index(1))).(featNamesTr{ifN}) = 1/sum([sTarget(:).nsamples])*...
                            (sTarget(1).nsamples*sTarget(1).(featNamesTr{ifN}) + sTarget(2).nsamples*sTarget(2).(featNamesTr{ifN}));
                        if any(strcmpi(fnTn,featNamesTr_SD{ifN}))
                            Targetnew(sindex(index(1))).(featNamesTr_SD{ifN}) = 1/sum([sTarget(:).nsamples])*...
                                (sTarget(1).nsamples*(sTarget(1).(featNamesTr_SD{ifN})^2+sTarget(1).(featNamesTr)^2)+...
                                sTarget(2).nsamples*(sTarget(2).(featNamesTr_SD{ifN})^2+sTarget(2).(featNamesTr{ifN})^2))-...
                                Targetnew(sindex(index(1))).(featNamesTr{ifN})^2;
                            Targetnew(sindex(index(1))).(featNamesTr_SD{ifN}) = sqrt(Targetnew(sindex(index(1))).(featNamesTr_SD{ifN}));
                        else
                            warning('no SD info')
                        end
                    else
                        warning('not all fieldnames present')
                    end
                end
                Targetnew(sindex(index(1))).('nsamples') = sTarget(1).nsamples+sTarget(2).nsamples;
            end
        end
        Targetnew(indextocut) = [];
        Target_tr = Targetnew;
        clearvars('Targetnew','sTarget','index','multiple','n','bin','sindex')

        %adjust weights as well
        switch fopt_objMethod
            case 'currenttraces'
                if contains(options_fopt.CT.method_CT,'weighted','IgnoreCase',true)
                    options_fopt.CT.wtr(indextocut,:) = [];
                end
            case 'features'
                options_fopt.FB.wtr(indextocut,:) = [];
            otherwise; error('why?');
        end
    end
end

if sort_flag
    [~,idx_s] = sort([Target_sp(:).Il]','descend');
    Target_sp = Target_sp(idx_s);
    Il_valsp = unique([Target_sp(:).Il],'stable');
    Vm_valsp = unique([Target_sp(:).Vm],'stable');

    %adjust weights as well
    if strcmpi(wGODA.type,'input')
        wGODA.w = wGODA.w(idx_s,:);
    end
    switch fopt_objMethod
        case 'currenttraces'
            if contains(options_fopt.CT.method_CT,'weighted','IgnoreCase',true)
                options_fopt.CT.w = options_fopt.CT.w(idx_s,:);
            end
        case 'features'
            options_fopt.FB.w = options_fopt.FB.w(idx_s,:);
        otherwise; error('why?');
    end

    if TauRecov_flag
        [~,idx_s] = sort([Target_tr(:).Il]','descend');
        Target_tr = Target_tr(idx_s);
        Il_valtr = unique([Target_tr(:).Il],'stable');
        Vm_valtr = unique([Target_tr(:).Vm],'stable');

        switch fopt_objMethod
            case 'currenttraces'
                if contains(options_fopt.CT.method_CT,'weighted','IgnoreCase',true)
                    options_fopt.CT.wtr = options_fopt.CT.wtr(idx_s,:);
                end
            case 'features'
                options_fopt.FB.wtr = options_fopt.FB.wtr(idx_s,:);
            otherwise; error('why?');
        end
    end
end


if plotFeat_flag
    plotFeatures(Target_sp,Il_valsp,Vm_valsp,featNamesT_plot,featNamesI_plot,Colors,errorbar_flag,...
        Target_tr,Il_valtr,Vm_valtr,plotTauRecov_flag,'figposition',fig_pos);
end


switch lower(IVdepend_targets)
    case 'features'
        [xdata,ydata,IlTau_val,options_GODA] = genXYFB(Target_sp,Vm_valsp,Il_valsp,Target_tr,Vm_valtr,Il_valtr,wGODA,options_GODA,...
            powera,powerb,recal_opt_FBIVdepend,TauRecov_flag,inclIratio_flag,corrTauDAIl0_flag);
        OinfDAinf_flag = 0;
    case 'currenttraces'

        error('not encoded')
        OinfDAinf_flag = 1;
    case 'postcurrenttraces'
        error('not encoded')
        data = load(['./Target/',filename]);
        OinfDAinf_flag = 1;
    otherwise
        error('invalid input fit Method')
end

if plotTau_flag
    plotTau(xdata,ydata,Vm_valsp,IlTau_val,errorbar_flag,'figposition',fig_pos)
end

% load fitDatabase
[fun,OinfIfun_names,TauOIfun_names,DAinfIfun_names,...
    TauDAIfun_names,TauOVfun_names,TauDAVfun_names,Grectfun_names, combTauO,combTauDA] = ...
    load_functionDB(varargin);


switch fitOrder
    case 'tO-tDA-GODA'
        % Fit TauO
        if display_flag
            fprintf('\n start fitting TauO \n')
        end
        if time_flag
            Timer(2).value = toc(starttimer);
            Timer(2).comment = 'startfitTauO';
        end

        [OutTauO,funO_Best,idxO_Best,fvalO_Best] = fitTau(xdata.TauO,ydata.TauO,'TauO',...
            fun,TauOIfun_names,TauOVfun_names,combTauO,options_TauO,plot_flag,IlTau_val,Vm_valsp,...
            Colors,type,memorySave_flag,plotTau_flag);
        % add to real output
        Out.TauO = OutTauO.TauO;
        fun_Best.TauO = funO_Best.TauO;
        idx_Best.TauO = idxO_Best.TauO;
        fval_Best.TauO = fvalO_Best.TauO;

        % Fit TauDA
        if display_flag
            fprintf('\n start fitting TauDA \n')
        end
        if time_flag
            Timer(3).value = toc(starttimer);
            Timer(3).comment = 'startfitTauDA';
        end
        [OutTauDA,funDA_Best,idxDA_Best,fvalDA_Best] = fitTau(xdata.TauDA,ydata.TauDA,'TauDA',...
            fun,TauDAIfun_names,TauDAVfun_names,combTauDA,options_TauDA,plot_flag,IlTau_val,Vm_valsp,...
            Colors,type,memorySave_flag,plotTau_flag);
        % add to real output
        Out.TauDA = OutTauDA.TauDA;
        fun_Best.TauDA = funDA_Best.TauDA;
        idx_Best.TauDA = idxDA_Best.TauDA;
        fval_Best.TauDA = fvalDA_Best.TauDA;

        % Fit Oinf DAinf Grect
        if display_flag
            fprintf('\n start fitting current dependecies Oinf, DAinf, Grect \n')
        end
        if time_flag
            Timer(4).value = toc(starttimer);
            Timer(4).comment = 'startfitI';
        end

        if OinfDAinf_flag
            error('not encoded')
        else

            [OutI,funI_Best,idxI_Best,fvalI_Best] = fitGODA_IpIssBased(xdata.Ipeak,ydata.Ipeak,xdata.Iss,ydata.Iss,xdata.Iratio,ydata.Iratio,...
                powera,powerb,fun_Best.TauO,fun_Best.TauDA,fun,OinfIfun_names,DAinfIfun_names,...
                Grectfun_names,options_GODA,plot_flag,display_flag,Il_valsp,Vm_valsp,Colors,type,inclIratio_flag,memorySave_flag,'figposition',fig_pos);
            % add to real output
            Out.I = OutI.I;
            fun_Best.Oinf = funI_Best.Oinf;
            fun_Best.DAinf = funI_Best.DAinf;
            fun_Best.Grect = funI_Best.Grect;
            idx_Best.I = idxI_Best.I;
            fval_Best.I = fvalI_Best.I;

        end
    case 'GODA-tDA-tO'

        % Fit Oinf DAinf Grect
        if display_flag
            fprintf('\n start fitting current dependecies Oinf, DAinf, Grect \n')
        end
        if time_flag
            Timer(2).value = toc(starttimer);
            Timer(2).comment = 'startfitI';
        end

        if OinfDAinf_flag
            error('not encoded')
        elseif inclIratio_flag
            error('incl Iratio not incoded for this fit order')
        else

            [OutI,funI_Best,idxI_Best,fvalI_Best] = fitGODA_allfeatBased(xdata.Ipeak,ydata.Ipeak,xdata.Iss,ydata.Iss,...
                ydata.TauOn,ydata.TauInact,powera,powerb,fun,OinfIfun_names,DAinfIfun_names,...
                Grectfun_names,options_GODA,plot_flag,display_flag,Il_valsp,Vm_valsp,Colors,type,memorySave_flag,'figposition',fig_pos);
            % add to real output
            Out.I = OutI.I;
            fun_Best.Oinf = funI_Best.Oinf;
            fun_Best.DAinf = funI_Best.DAinf;
            fun_Best.Grect = funI_Best.Grect;
            idx_Best.I = idxI_Best.I;
            fval_Best.I = fvalI_Best.I;
        end

        % Fit TauDA
        if display_flag
            fprintf('\n start fitting TauDA \n')
        end
        if time_flag
            Timer(3).value = toc(starttimer);
            Timer(3).comment = 'startfitTauDA';
        end
        [OutTauDA,funDA_Best,idxDA_Best,fvalDA_Best] = fitTau(xdata.TauDA,ydata.TauDA,'TauDA',...
            fun,TauDAIfun_names,TauDAVfun_names,combTauDA,options_TauDA,plot_flag,IlTau_val,Vm_valsp,...
            Colors,type,memorySave_flag,plotTau_flag);
        % add to real output
        Out.TauDA = OutTauDA.TauDA;
        fun_Best.TauDA = funDA_Best.TauDA;
        idx_Best.TauDA = idxDA_Best.TauDA;
        fval_Best.TauDA = fvalDA_Best.TauDA;

         % Fit TauO
        if display_flag
            fprintf('\n start fitting TauO \n')
        end
        if time_flag
            Timer(4).value = toc(starttimer);
            Timer(4).comment = 'startfitTauO';
        end

        [OutTauO,funO_Best,idxO_Best,fvalO_Best] = fitTau(xdata.TauO,ydata.TauO,'TauO',...
            fun,TauOIfun_names,TauOVfun_names,combTauO,options_TauO,plot_flag,IlTau_val,Vm_valsp,...
            Colors,type,memorySave_flag,plotTau_flag,...
            'nonlcon',1,'TauDA',fun_Best.TauDA,'DAinf',fun_Best.TauDA,'powera',powera,'powerb',powerb);
        % add to real output
        Out.TauO = OutTauO.TauO;
        fun_Best.TauO = funO_Best.TauO;
        idx_Best.TauO = idxO_Best.TauO;
        fval_Best.TauO = fvalO_Best.TauO;
    case 'noFirstFit'
    otherwise
        error('wrong fit order')
end
% run ODE
%warning('run ode off always')
if plot_flag
    %time consuming
    runODE(Vm_valsp,Target_sp,OdeOpts,fun_Best,powera,powerb,Colors,plotTauRecov_flag,Vm_valtr,Target_tr)
end


% Final optimization
if display_flag
    fprintf('\n final optimization \n')
end
if time_flag
Timer(5).value = toc(starttimer);
Timer(5).comment = 'optimize all';
end

[Out_final,funfopt_Best,idxfopt_Best,fvalfopt_Best] = final_opt(fopt_select,fopt_objMethod,fopt_solverMethod,...
    OinfIfun_names,TauOIfun_names,DAinfIfun_names,TauDAIfun_names,TauOVfun_names,TauDAVfun_names,Grectfun_names,combTauO,combTauDA,...
    Out,fun,powera,powerb,Vm_valsp,Il_valsp,idx_Best,xdata,ydata,Target_sp,Target_tr,Vm_valtr,options_fopt,...
    display_flag,plot_flag,plotTauRecov_flag,TauRecov_flag,Colors,memorySave_flag);
Out.(fopt_select) = Out_final.(fopt_select);

if plot_flag
    runODE(Vm_valsp,Target_sp,OdeOpts,funfopt_Best,powera,powerb,Colors,plotTauRecov_flag,Vm_valtr,Target_tr)
end
if time_flag
Timer(6).value = toc(starttimer);
Timer(6).comment = 'end';
Out.Timer = Timer;
end
end