function runODE(Vm_val,Target,OdeOpts_set,fun_Best,powera,powerb,Colors,plotTauRecov_flag,Vm_valtr,Target_tr,varargin)
fignr = get(gcf,'number')+1;

if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'fignr'))
            fignr = varargin{find(strcmpi(varargin,'fignr'))+1};
        end
    end
end

for iv=1:length(Vm_val)
    idxv = find([Target(:).Vm]==Vm_val(iv));
    Il_valv = [Target(idxv).Il];
    for iIl = 1:length(Il_valv)

        idxt = idxv(iIl);
        OSpstart=Target(idxt).OSpstart; OSpd=Target(idxt).OSpd;OSipa=Target(idxt).Il;
        if isfield(Target,'time')
            tend = Target(idxt).time(end);
        else
            tend = 2*(OSpstart+OSpd);
        end
        OSstep = @(t) double(t>=OSpstart&t<=OSpd+OSpstart);
        OSi = @(t) OSipa*OSstep(t);  % [W/m^2]
        OdeOpts = odeset(OdeOpts_set{:});
        [t,m]=ode15s(@(t,y) States_Vclamp(t,y(1),y(2),fun_Best.Oinf,fun_Best.TauO,...
            fun_Best.DAinf,fun_Best.TauDA,OSi,Target(idxt).Vm),[0,tend],[0,1],OdeOpts);
        iChR2=fun_Best.Grect([Target(idxt).Vm.*ones(size(t)),OSi(t)]).*m(:,1).^powera.*(m(:,2)).^powerb;

        figure(fignr)
        subplot(ceil(length(Vm_val)/2),2,iv)
        plot(t,iChR2,'color',Colors(iIl),'DisplayName',['Fit: ',num2str(Target(idxt).Il)]);
        hold on
        try
        plot(Target(idxt).time,Target(idxt).I,'--','color',Colors(iIl),'DisplayName',['raw data']);
        end
        try
        plot(OSpstart+OSpd,Target(idxt).Iss,'*','color',Colors(iIl))
        end
        try
        plot(OSpstart:0.01:OSpstart+0.1,ones(size(OSpstart:0.01:OSpstart+0.1)).*Target(idxt).Ipeak,'.','color',Colors(iIl))
        end
    end
    hold off
    title(num2str(Target(idxt).Vm))
end

if plotTauRecov_flag
    for iv=1:length(Vm_valtr)
    idxv = find([Target_tr(:).Vm]==Vm_valtr(iv));
    Il_valv = [Target_tr(idxv).Il];
    for iIl = 1:length(Il_valv)

        idxt = idxv(iIl);
        OSpstart=Target_tr(idxt).Pstart; PD=Target_tr(idxt).PD;OSipa=Target_tr(idxt).Il;
        intervals = Target_tr(idxt).Intervals;
        for iInt = 1:length(intervals)
            OSprp = intervals(iInt)+PD;      % Pulse repetition period (s)
            OSdc = PD/OSprp;
            OSpd= 2*OSprp;
            OSstep = @(t) double(mod(t-OSpstart,OSprp)<=OSdc*OSprp).*double(t>=OSpstart&t<=OSpd+OSpstart);
            if isfield(Target_tr(idxt),'time')
                tend_tr = Target_tr(idxt).time{iInt}(end);
            else
                tend_tr = OSpstart*2+OSpd;
            end
            OSi = @(t) OSipa*OSstep(t);  % [W/m^2]
            OdeOpts = odeset(OdeOpts_set{:});
            [t,m]=ode15s(@(t,y) States_Vclamp(t,y(1),y(2),fun_Best.Oinf,fun_Best.TauO,...
                fun_Best.DAinf,fun_Best.TauDA,OSi,Target_tr(idxt).Vm),[0,tend_tr],[0,1],OdeOpts);
            iChR2=fun_Best.Grect([Target_tr(idxt).Vm.*ones(size(t)),OSi(t)]).*m(:,1).^powera.*(m(:,2)).^powerb;

            figure(fignr+1)
            subplot(ceil(length(Vm_valtr)/2),2,iv)
            plot(t,iChR2,'color',Colors(iIl),'DisplayName',['Fit: ',num2str(Target_tr(idxt).Il)]);
            hold on
            try
                plot(Target_tr(idxt).time{iInt},Target_tr(idxt).I{iInt},'--','color',Colors(iIl),'DisplayName',['raw data']);
            end
        end
    end
    hold off
    title(num2str(Target_tr(idxt).Vm))
    end
end
end