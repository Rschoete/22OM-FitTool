function [xdata,ydata,Il_valTAU,options_GODA] = genXYFB(Target_sp,Vm_valsp,Il_valsp,Target_tr,...
    Vm_valtr,Il_valtr,wGODA,options_GODA,powera,powerb,recal_opt,TauRecov_flag,inclIratio_flag,corrTauDAIl0_flag)

fnTsp = fieldnames(Target_sp);
if TauRecov_flag
fnTtr = fieldnames(Target_tr);
end

% Idata
xdata.Ipeak = [[Target_sp(:).Vm]',[Target_sp(:).Il]'];
xdata.Iss = [[Target_sp(:).Vm]',[Target_sp(:).Il]'];
ydata.Ipeak = [Target_sp(:).Ipeak]';
ydata.Iss = [Target_sp(:).Iss]';
idx = isnan(ydata.Ipeak)|isnan(ydata.Iss);
xdata.Ipeak(idx,:) = [];
ydata.Ipeak(idx) = [];
xdata.Iss(idx,:) = [];
ydata.Iss(idx) = [];

xdata.Iratio = [[Target_sp(:).Vm]',[Target_sp(:).Il]'];
ydata.Iratio = [Target_sp(:).Iratio]';
idxr = isnan(ydata.Iratio);
xdata.Iratio(idxr,:) = [];
ydata.Iratio(idxr) = [];

if any(strcmpi(fnTsp,'Ipeak_SD'))
    ydata.Ipeak_SD = [Target_sp(:).Ipeak_SD]';
    ydata.Ipeak_SD(idx) = [];
end
if any(strcmpi(fnTsp,'Iss_SD'))
    ydata.Iss_SD = [Target_sp(:).Iss_SD]';
    ydata.Iss_SD(idx) = [];
end
if any(strcmpi(fnTsp,'Iratio_SD'))
    ydata.Iratio_SD = [Target_sp(:).Iratio_SD]';
    ydata.Iratio_SD(idxr) = [];
end
% create correct weight arrays
switch lower(wGODA.type)
    case 'sd'
        options_GODA.wpeak = 1./[ydata.Ipeak_SD];
        options_GODA.wss = 1./[ydata.Iss_SD];
        if inclIratio_flag
            options_GODA.wratio = 1./[ydata.Iratio_SD];
        end
    case 'rel'
        options_GODA.wpeak = 1./[ydata.Ipeak];
        options_GODA.wss = 1./[ydata.Iss];
        if inclIratio_flag
            options_GODA.wratio = 1./[ydata.Iratio];
        end
    case 'input'
        options_GODA.wpeak = wGODA.w(:,1); options_GODA.wpeak(idx) = [];
        options_GODA.wss = wGODA.w(:,2); options_GODA.wss(idx) = [];
        if inclIratio_flag
            options_GODA.wratio = wGODA.w(:,3); options_GODA.wratio(idxr) = [];
        end
    case 'equal'
        options_GODA.wpeak = ones(size(ydata.Ipeak,1),1);
        options_GODA.wss = ones(size(ydata.Iss,1),1);
        if inclIratio_flag
            options_GODA.wratio = ones(size(ydata.Iratio,1),1);
        end
    otherwise
        error('false input')
end

% TauO
xdata.TauO=[[Target_sp(:).Vm]',[Target_sp(:).Il]'];
xdata.TauOn = xdata.TauO;
for i=1:length(Vm_valsp)
    idx_Vm = [Target_sp.Vm]==Vm_valsp(i);
    AvgTauOff(i)=1/sum([Target_sp(idx_Vm).nsamples])*...
        sum([Target_sp(idx_Vm).nsamples].*[Target_sp(idx_Vm).TauOff]);

    if any(strcmpi(fnTsp,'TauOff_SD'))

        AvgTauOff_SD(i) = 1/sum([Target_sp(idx_Vm).nsamples])*...
            sum([Target_sp(idx_Vm).nsamples].*...
            ([Target_sp(idx_Vm).TauOff_SD].^2+...
            [Target_sp(idx_Vm).TauOff].^2))-...
            AvgTauOff(i)^2;
        AvgTauOff_SD(i) = sqrt(AvgTauOff_SD(i));
    end
end
xdata.TauO=[xdata.TauO;[Vm_valsp',zeros(size(Vm_valsp))']];
xdata.TauOff = [Vm_valsp',zeros(size(Vm_valsp))'];
ydata.TauO=[[Target_sp(:).TauOn]';AvgTauOff'];
idx = isnan(ydata.TauO);
xdata.TauO(idx,:) = [];
ydata.TauO(idx) = [];

ydata.TauOn = [Target_sp(:).TauOn]';
ydata.TauOff = AvgTauOff';
idx_tOn = isnan(ydata.TauOn);
xdata.TauOn(idx_tOn,:) = [];
ydata.TauOn(idx_tOn,:) = [];
idx_tOff = isnan(ydata.TauOff);
xdata.TauOff(idx_tOff,:) = [];
ydata.TauOff(idx_tOff,:) = [];


if any(strcmpi(fnTsp,'TauOn_SD'))
    ydata.TauO_SD = [[Target_sp(:).TauOn_SD]'];

else
    ydata.TauO_SD = [];
end

if any(strcmpi(fnTsp,'TauOff_SD'))
    ydata.TauO_SD = [ydata.TauO_SD;AvgTauOff_SD'];
    ydata.TauO_SD(idx) = [];

end
if powera~=1
    ydata.TauO=recalTauO(ydata.TauO,powera,xdata.TauO,recal_opt.TauO);
end

%TauDA
xdata.TauDA=[[Target_sp(:).Vm]',[Target_sp(:).Il]'];
xdata.TauInact = [[Target_sp(:).Vm]',[Target_sp(:).Il]'];
if TauRecov_flag
    for i=1:length(Vm_valtr)
        idx_Vm = [Target_tr.Vm]==Vm_valtr(i);
        AvgTauRecov(i) = 1/sum([Target_tr(idx_Vm).nsamples])*...
            sum([Target_tr(idx_Vm).nsamples].*[Target_tr(idx_Vm).TauRecov]);
        if any(strcmpi(fnTtr,'TauRecov_SD'))
            AvgTauRecov_SD(i) = 1/sum([Target_tr(idx_Vm).nsamples])*...
                sum([Target_tr(idx_Vm).nsamples].*...
                ([Target_tr(idx_Vm).TauRecov_SD].^2+...
                [Target_tr(idx_Vm).TauRecov].^2))-...
                AvgTauOff(i)^2;
            AvgTauRecov_SD(i) = sqrt(AvgTauRecov_SD(i));

        end
    end
    if corrTauDAIl0_flag
        % correction for TauDAIlO not exactly equal to Tau recov
        Vm_valIr = unique(xdata.Iratio(:,1),'stable');
        [iVm_valIrTr] = intersect(Vm_valtr,Vm_valIr,'stable'); %orderd according to appearance in Vm_valtr
        if numel(iVm_valIrTr) == numel(Vm_valtr)
            for iiVm = 1:length(iVm_valIrTr)
                idx_VmTr = Vm_valtr==iVm_valIrTr(iiVm);
                idx_VmIr = xdata.Iratio(:,1)==iVm_valIrTr(iiVm);
                CF = 1-log(1/(1-median(ydata.Iratio(idx_VmIr))));
                if CF<0.1
                    CF = 0.1;
                    warning('CF for Tau Recov fixed to 0.1')
                end
                AvgTauRecov(idx_VmTr) = AvgTauRecov(idx_VmTr)/CF;
                if any(strcmpi(fnTtr,'TauRecov_SD'))
                    AvgTauRecov_SD(idx_VmTr) = AvgTauRecov_SD(idx_VmTr)/CF;
                end
            end
        else
            CF = 1-log(1/(1-median(ydata.Iratio)));
            if CF<0.1
                CF = 0.1;
                warning('CF for Tau Recov fixed to 0.1')
            end
            AvgTauRecov = AvgTauRecov/CF;
            if any(strcmpi(fnTtr,'TauRecov_SD'))
                AvgTauRecov_SD = AvgTauRecov_SD/CF;
            end

        end

    end
    xdata.TauDA=[xdata.TauDA;[Vm_valtr',zeros(size(Vm_valtr))']];
    xdata.TauRecov = [Vm_valtr',zeros(size(Vm_valtr))'];
    ydata.TauDA=[[Target_sp(:).TauInact]';AvgTauRecov'];
    idx = isnan(ydata.TauDA);
    xdata.TauDA(idx,:) = [];
    ydata.TauDA(idx) = [];

    ydata.TauInact = [Target_sp(:).TauInact]';
    ydata.TauRecov = AvgTauRecov';

    idx_ti = isnan(ydata.TauInact);
    xdata.TauInact(idx_ti,:) = [];
    ydata.TauInact(idx_ti,:) = [];

    idx_tr = isnan(ydata.TauRecov);
    xdata.TauInact(idx_tr,:) = [];
    ydata.TauInact(idx_tr,:) = [];

    if any(strcmpi(fnTsp,'TauInact_SD'))
        ydata.TauDA_SD=[[Target_sp(:).TauInact_SD]'];
    else
        ydata.TauDA_SD=[];
    end

    if any(strcmpi(fnTtr,'TauRecov_SD'))
        ydata.TauDA_SD=[ydata.TauDA_SD;AvgTauRecov_SD'];
        ydata.TauDA_SD(idx) = [];
    end
else

    ydata.TauDA=[Target_sp(:).TauInact]';
    idx = isnan(ydata.TauDA);
    xdata.TauDA(idx,:) = [];
    ydata.TauDA(idx) = [];

    ydata.TauInact = [Target_sp(:).TauInact]';

    idx_ti = isnan(ydata.TauInact);
    xdata.TauInact(idx_ti,:) = [];
    ydata.TauInact(idx_ti,:) = [];

    if any(strcmpi(fnTsp,'TauInact_SD'))
        ydata.TauDA_SD=[Target_sp(:).TauInact_SD]';
        ydata.TauDA_SD(idx) = [];
    end
    if corrTauDAIl0_flag
        warning('corrTauDAIl0_flag is on but Tau recov not included == useless')
    end
end
if powerb ~= 1
    ydata.TauI=recalTauDA(ydata.TauDA,powerb,xdata.TauDA,recal_opt.TauDA);%[Target(:).Iss],[Target(:).Ipeak]);
end
Il_valTAU =[Il_valsp,0];

% when higher power orders are requested, the targets need to be changed as
% the targets are based on a monoexponential fit
    function ym = recalTauO(tTau,powera,tX,options)

        ym=10*ones(size(tTau));
        for io=1:size(tTau,1)
            if tX(io,2)~=0
                X0 = tTau(io)/((powera-1)/2+1);
                if isfield(options,'dt') && isfield(options,'tend')
                    to=0:options.dt:options.tend;
                else
                    dt = X0/100; tend = X0*20;
                    to = 0:dt:tend;
                end
                if isfield(options,'ub')
                    UB = options.ub;
                else
                    UB = X0*20;
                end
                if isfield(options,'lb')
                    LB = options.lb;
                else
                    LB = 0;
                end
                ty=(1-exp(-to./tTau(io)));
                OBJo=@(p2,t) (1-exp(-t/p2)).^powera;
                opt=optimoptions(options.optimoptions{:});
                problemo=createOptimProblem('lsqcurvefit',...
                    'objective',OBJo,...
                    'xdata',to,'ydata',ty,...
                    'x0',X0,...
                    'lb',LB,...
                    'ub',UB, 'options',opt);
                mso = MultiStart(options.msopt{:});
                [ym(io)]= run(mso,problemo,options.msnr);
            end
        end
        % TauOff is different
        idxI0=tX(:,2)==0;
        ym(idxI0)=tTau(idxI0)*powera;
    end

    function ym = recalTauDA(tTau,powerb,tX,options)

        ym=10*ones(size(tTau));
        for io=1:size(tTau,1)
            X0 = tTau(io)/((powerb-1)/2+1);
            if tX(io,2)~=0
                if isfield(options,'dt') && isfield(options,'tend')
                    to=0:options.dt:options.tend;
                else
                    to = 0:X0/100:X0*20;
                end
                if isfield(options,'ub')
                    UBo = options.ub;
                else
                    UBo = X0*20;
                end
                if isfield(options,'lb')
                    LBo = options.lb;
                else
                    LBo = 0;
                end
            else
                if isfield(options,'dtr') && isfield(options,'tendr')
                    to = 0:options.dtr:options.tendr;
                else
                    to = 0:X0/100:X0*20;
                end
                if isfield(options,'ubr')
                    UBo = options.ubr;
                else
                    UBo = X0*20;
                end
                if isfield(options,'lbr')
                    LBo = options.lbr;
                else
                    LBo = 0;
                end
            end
                ty=(1-exp(-to./tTau(io)));
                OBJo=@(p2,t) (1-exp(-t/p2)).^powerb;
                opt=optimoptions(options.optimoptions{:});
                problemo=createOptimProblem('lsqcurvefit',...
                    'objective',OBJo,...
                    'xdata',to,'ydata',ty,...
                    'x0',X0,...
                    'lb',LBo,...
                    'ub',UBo, 'options',opt);
                mso = MultiStart(options.msopt{:});
                [ym(io)]= run(mso,problemo,options.msnr);
        end
    end
end