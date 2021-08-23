function [Out_iChR2_sp,time_iChR2_sp,Out_iChR2_tp,time_iChR2_tp,Ipeak,Iss,Iratio,TauOn,TauInact,TauOff,TauRecov] = ...
    collDATA_fittool(Input,Target,intm_flag,tp_flag,OdeSet,nrruns,PLOT,varargin)

Target_sp = Target.SingleP;
Target_tr = Target.TauRecov;

for ifit = 1:1+intm_flag
    if ifit==1
        powera = Input.all.powera;
        powerb = Input.all.powerb;

        Oinf = str2func(Input.all.param.OIfun);
        Oinf = @(X) Oinf(Input.all.param.pOI,X(:,2));
        DAinf = str2func(Input.all.param.DAIfun);
        DAinf = @(X) DAinf(Input.all.param.pDAI,X(:,2));

        combTauO = str2func(Input.all.param.combTauO);
        TauOI = str2func(Input.all.param.TauOIfun);
        TauOI = @(X) TauOI(Input.all.param.pTOI,X(:,2));
        TauOV = str2func(Input.all.param.TauOVfun);
        TauOV = @(X) TauOV(Input.all.param.pTOV,X(:,1));
        TauO = @(X) combTauO(TauOV(X),TauOI(X));

        combTauDA = str2func(Input.all.param.combTauDA);
        TauDAI = str2func(Input.all.param.TauDAIfun);
        TauDAI = @(X) TauDAI(Input.all.param.pTDAI,X(:,2));
        TauDAV = str2func(Input.all.param.TauDAVfun);
        TauDAV = @(X) TauDAV(Input.all.param.pTDAV,X(:,1));
        TauDA = @(X) combTauDA(TauDAV(X),TauDAI(X));

        Grect = str2func(Input.all.param.Gfun);
        Grect = @(X) Grect(Input.all.param.pG,X(:,1));
    else
        powera = Input.I.param.powera;
        powerb = Input.I.param.powerb;

        Oinf = str2func(Input.I.param.OIfun);
        Oinf = @(X) Oinf(Input.I.param.pOI,X(:,2));
        DAinf = str2func(Input.I.param.DAIfun);
        DAinf = @(X) DAinf(Input.I.param.pDAI,X(:,2));

        combTauO = str2func(Input.TauO.param.combTau);
        TauOI = str2func(Input.TauO.param.TauIfun);
        TauOI = @(X) TauOI(Input.TauO.param.pI,X(:,2));
        TauOV = str2func(Input.TauO.param.TauVfun);
        TauOV = @(X) TauOV(Input.TauO.param.pV,X(:,1));
        TauO = @(X) combTauO(TauOV(X),TauOI(X));

        combTauDA = str2func(Input.TauDA.param.combTau);
        TauDAI = str2func(Input.TauDA.param.TauIfun);
        TauDAI = @(X) TauDAI(Input.TauDA.param.pI,X(:,2));
        TauDAV = str2func(Input.TauDA.param.TauVfun);
        TauDAV = @(X) TauDAV(Input.TauDA.param.pV,X(:,1));
        TauDA = @(X) combTauDA(TauDAV(X),TauDAI(X));

        Grect = str2func(Input.I.param.Gfun);
        Grect = @(X) Grect(Input.I.param.pG,X(:,1));
    end

for iTsp = 1:length(Target_sp)

    Vm = Target_sp(iTsp).Vm;
    OSpstart = Target_sp(iTsp).OSpstart; OSpd = Target_sp(iTsp).OSpd; OSipa = Target_sp(iTsp).Il;
    if isfield(Target_sp,'time')
        tend = Target_sp(iTsp).time(end);
    else
        tend = 2*(OSpstart+OSpd);
    end
    OSstep = @(t) double(t>=OSpstart&t<=OSpd+OSpstart);
    OSi = @(t) OSipa*OSstep(t);  % [W/m^2]
    OdeOpts = odeset(OdeSet{:});
    [t,m]=ode15s(@(t,y) States_Vclamp(t,y(1),y(2),Oinf,TauO,...
        DAinf,TauDA,OSi,Vm),[0,tend],[0,1],OdeOpts);
    iChR2=Grect([Vm.*ones(size(t)),OSi(t)]).*m(:,1).^powera.*(m(:,2)).^powerb;

    features = Extract_feat(t,iChR2,OSpd,OSpstart,nrruns,PLOT,varargin);
    % output is downsampled to reduce memory  usage
    downsf = max(floor(length(iChR2)/3e3),1);
    Out_iChR2_sp{ifit,iTsp} = iChR2(1:downsf:end);
    time_iChR2_sp{ifit,iTsp} = t(1:downsf:end);
    Ipeak(ifit,iTsp) = features.Ipeak;
    Iss(ifit,iTsp) = features.Iss;
    Iratio(ifit,iTsp) = features.Iratio;
    TauOn(ifit,iTsp) = features.TauOn;
    TauInact(ifit,iTsp) = features.TauInact;
    TauOff(ifit,iTsp) = features.TauOff;
end

if tp_flag

    for iTtp = 1:length(Target_tr)
        intervals = Target_tr(iTtp).Intervals;
        Peak1 = zeros(length(intervals),1);
        Peak2 = zeros(length(intervals),1);
        Ratio = zeros(length(intervals),1);
        for iInt = 1:length(intervals)
            Vm = Target_tr(iTtp).Vm;
            OSpstart = Target_tr(iTtp).Pstart; OSpd = Target_tr(iTtp).PD; OSipa = Target_tr(iTtp).Il;
            OSprp = OSpd + intervals(iInt);
            OSdc = OSpd/OSprp;
            if isfield(Target_tr,'tend')
                tend = Target_tr(iTtp).tend(iInt);
            else
                tend = 2*(OSpstart+OSprp);
            end
            OSstep = @(t) double(mod(t-OSpstart,OSprp)<=OSdc*OSprp).*double(t>=OSpstart&t<=2*OSprp);
            OSi = @(t) OSipa*OSstep(t);  % [W/m^2]
            OdeOpts = odeset(OdeSet{:});
            [t,m]=ode15s(@(t,y) States_Vclamp(t,y(1),y(2),Oinf,TauO,...
                DAinf,TauDA,OSi,Vm),[0,tend],[0,1],OdeOpts);
            iChR2=Grect([Vm.*ones(size(t)),OSi(t)]).*m(:,1).^powera.*(m(:,2)).^powerb;
            downsf = max(floor(length(iChR2)/1e5),1);
            Out_iChR2_tp{ifit,iTtp,iInt} = iChR2(1:downsf:end);
            time_iChR2_tp{ifit,iTtp,iInt} = t(1:downsf:end);

            loc1=t>=OSpstart & t<=(OSpstart+OSpd+1/2*intervals(iInt));
            loc2=t>=OSpstart+OSprp & t<=OSpstart+OSprp+OSpd+1/2*intervals(iInt);
            Peak1(iInt)=max(abs(iChR2(loc1)));
            %Peak2(i)=max(abs(Out(i).iChR2(loc2)));
            try
                Peaks2 = findpeaks(abs(iChR2(loc2)));
                Peak2(iInt) = Peaks2(1);
            catch
                Peak2(iInt) = 0;
            end

            Ratio(iInt)=Peak2(iInt)/Peak1(iInt);

        end

        xx = [0:1e-3:intervals(end)];
        RATIO=spline(intervals,Ratio,xx);
        [~,idx]=min(abs(RATIO-(1-exp(-1))));
        if any(Ratio<=(1-exp(-1))) && any(Ratio>=(1-exp(-1)))
            TauRecov(ifit,iTtp)=xx(idx(1));
        else
            TauRecov(ifit,iTtp)=nan;
        end
    end
end
end
end