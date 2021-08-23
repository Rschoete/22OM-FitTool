function [Out,fun_Best,idx_Best,fval_Best] = fitGODA_allfeatBased(xIp,yIp,xIss,yIss,yton,ytinact,powera,powerb,fun,OinfIfun_names,DAinfIfun_names,...
    Grectfun_names,options,plot_flag,display_flag,Il_valsp,Vm_valsp,Colors,type,memorySave_flag,varargin)

if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'figposition'))
            figposition = varargin{find(strcmpi(varargin,'figposition'))+1};
            figpos_flag = 1;
        end
    end
end

if any(xIp(:)-xIss(:))
    error('should be equal')
end
if length(OinfIfun_names)~=length(DAinfIfun_names) || length(OinfIfun_names)~=length(Grectfun_names)
    error('fun_names should have same size')
end

fignr = get(gcf,'number')+1;
ifn = 0;

% make sure targets are in column format
if ~iscolumn(yIp)
    if isvector(yIp)
        yIp = yIp';
    else
        error('why yIp  matrix?')
    end
end
if ~iscolumn(yIss)
    if isvector(yIss)
        yIss = yIss';
    else
        error('why yIp matrix?')
    end
end
if ~iscolumn(yton)
    if isvector(yton)
        yton = yton';
    else
        error('why yton matrix?')
    end
end
if ~iscolumn(ytinact)
    if isvector(ytinact)
        ytinact = ytinact';
    else
        error('why ytinact matrix?')
    end
end

if any([size(yton,1),size(ytinact,1),size(yIss,1)]-size(yIp,1))
    error('size of yton, ytinact, yIss , yIp has to be the same otherwise this method not possible')
end
for ifn=1:length(OinfIfun_names)

    if display_flag
        fprintf('\n %d out off %d \n',ifn,length(OinfIfun_names))
    end

    OinfIfun = fun.(OinfIfun_names{ifn}).fun;
    Pidx(1) = fun.(OinfIfun_names{ifn}).nrpara;
    LBOinfI = fun.(OinfIfun_names{ifn}).LB;
    UBOinfI = fun.(OinfIfun_names{ifn}).UB;
    X0OinfI = fun.(OinfIfun_names{ifn}).X0;
    VarOinfI = fun.(OinfIfun_names{ifn}).var;

    DAinfIfun = fun.(DAinfIfun_names{ifn}).fun;
    Pidx(2) = fun.(DAinfIfun_names{ifn}).nrpara+Pidx(1);
    LBDAinfI = fun.(DAinfIfun_names{ifn}).LB;
    UBDAinfI = fun.(DAinfIfun_names{ifn}).UB;
    X0DAinfI = fun.(DAinfIfun_names{ifn}).X0;
    VarDAinfI = fun.(DAinfIfun_names{ifn}).var;

    Grectfun = fun.(Grectfun_names{ifn}).fun;
    Pidx(3) = fun.(Grectfun_names{ifn}).nrpara+Pidx(2);
    LBGrect = fun.(Grectfun_names{ifn}).LB;
    UBGrect = fun.(Grectfun_names{ifn}).UB;
    X0Grect = fun.(Grectfun_names{ifn}).X0;
    VarGrect = fun.(Grectfun_names{ifn}).var;

    X0 = [X0OinfI,X0DAinfI,X0Grect];
    LB = [LBOinfI,LBDAinfI,LBGrect];
    UB = [UBOinfI,UBDAinfI,UBGrect];

    pOIidx = 1:Pidx(1);
    pDAIidx = Pidx(1)+1:Pidx(2);
    pGidx = Pidx(2)+1:Pidx(3);

    yfunx = @(t,taux,xinf,x0) xinf+(x0-xinf).*exp(-t./taux);
    TauO = @(X) yton(xIp(:,1)==X(:,1)&xIp(:,2)==X(:,2));
    TauDA = @(X) ytinact(xIss(:,1)==X(:,1)&xIss(:,2)==X(:,2));

    Oton = @(p,X,t) yfunx(t,TauO(X),OinfIfun(p(pOIidx),X(:,VarOinfI)),0);
    DAton = @(p,X,t) yfunx(t,TauDA(X),DAinfIfun(p(pDAIidx),X(:,VarDAinfI)),1);

    iChR2On = @(p,X,t) Grectfun(p(pGidx),X(:,VarGrect)).*(Oton(p,X,t)).^powera.*(DAton(p,X,t)).^powerb;

    OBJ=@(p) objfunIpeakIss(p,xIp,yIp,yIss,iChR2On,options);

    XO = options.I_nonlcon;
    XDA = options.I_nonlcon;
    XV1 = options.V1_nonlcon; if isrow(XV1); XV1 = XV1(:); end
    XV2 = options.V2_nonlcon; if isrow(XV2); XV2 = XV2(:); end
    XO0 = 0; if numel(VarOinfI)==2; XO0 = [[XV1;XV2],zeros(size([XV1;XV2]))]; end
    XDA0 = 0; if numel(VarDAinfI)==2; XDA0 = [[XV1;XV2],zeros(size([XV1;XV2]))]; end
    if numel(VarOinfI)==2 && size(XO,2)~=2
        XO = XO(:)';
        Vvals = [XV1(:);XV2(:)];
        [Vvals,XO] = meshgrid(Vvals,XO);
        XO = [Vvals(:),XO(:)];

    elseif numel(VarOinfI)==1 &&  size(XO,2)==2
        error('check input nonlincon both IV but funtion single column required')
    elseif numel(VarOinfI)==1 &&  isrow(XO)
        XO = XO(:);
    end
    if numel(VarDAinfI)==2 && size(XDA,2)~=2

        XDA = XDA(:)';
        Vvals = [XV1(:);XV2(:)];
        [Vvals,XDA] = meshgrid(Vvals,XDA);
        XDA = [Vvals(:),XDA(:)];
    elseif numel(VarDAinfI)==1 &&  size(XDA,2)==2
        error('check input nonlincon both IV but funtion single column required')
    elseif numel(VarDAinfI)==1 &&  isrow(XDA)
        XDA = XDA(:);
    end

    nonlcon = @(p) nonlincon(p,powera,powerb,OinfIfun,pOIidx,...
        DAinfIfun,pDAIidx,Grectfun,pGidx,...
        XO,XO0,XDA,XDA0,XV1,XV2,options.extrac_nonlcon,options.extraceq_nonlcon,...
        [],[],[]);

    switch lower(options.solverMethod)
        case 'ms'
            disp('start ms GODA')
            opt=optimoptions(options.optimoptions{:});
            problem=createOptimProblem('fmincon',...
                'objective',OBJ,...
                'nonlcon',nonlcon,...
                'x0',X0,...
                'lb',LB,...
                'ub',UB,...
                'options',opt);
            ms = MultiStart(options.msopt{:});
            [p,fval]= run(ms,problem,options.msnr);

        case 'psobc'
            disp('start psoBC GODA')
            problem = struct();
            problem = options.psoBC; % add simulation options
            if ~isfield(problem.options,'PopulationSize')
                problem.options.PopulationSize = length(LB)*10;
            end
            if isfield(options.psoBC.options,'PlotFcns')
                figure
            end
            problem.fitnessfcn = OBJ;
            problem.nvars = length(LB);
            problem.UB = UB;
            problem.LB = LB;
            problem.nonlcon = nonlcon;

            [p,fval] = pso(problem);

        case 'ga'
            disp('start ga GODA')
            problem = struct();
            problem = options.ga; % add simulation options
            if ~isfield(problem.options,'PopulationSize')
                problem.options.PopulationSize = length(LB)*10;
            end
            if isfield(options.ga.options,'PlotFcns')
                figure
            end
            problem.fitnessfcn = OBJ;
            problem.nvars = length(LB);
            problem.ub = UB;
            problem.lb = LB;
            problem.nonlcon = nonlcon;
            problem.solver = 'ga';

            [p,fval] = ga(problem);
    end

    Out.I.param(ifn).p = p;
    Out.I.param(ifn).pOI = p(pOIidx);
    Out.I.param(ifn).pOIidx = pOIidx;
    Out.I.param(ifn).pDAI = p(pDAIidx);
    Out.I.param(ifn).pDAIidx = pDAIidx;
    Out.I.param(ifn).pG = p(pGidx);
    Out.I.param(ifn).pGidx = pGidx;
    Out.I.param(ifn).fval = fval;
    Out.I.param(ifn).OIn = OinfIfun_names{ifn};
    Out.I.param(ifn).OIfun = func2str(OinfIfun);
    Out.I.param(ifn).DAIn = DAinfIfun_names{ifn};
    Out.I.param(ifn).DAIfun = func2str(DAinfIfun);
    Out.I.param(ifn).Gn = Grectfun_names{ifn};
    Out.I.param(ifn).Gfun = func2str(Grectfun);
    Out.I.param(ifn).powera = powera;
    Out.I.param(ifn).powerb = powerb;
    if ~memorySave_flag
        Out.I.param(ifn).X0 = X0;
        Out.I.param(ifn).LB = LB;
        Out.I.param(ifn).UB = UB;
        Out.I.xIp = xIp;
        Out.I.yIp = yIp;
        Out.I.xIss = xIss;
        Out.I.yIss = yIss;
        Out.I.options = options;
    end

    if ifn==1

        fval_Best.I = fval;
        idx_Best.I = [ifn];
        fun_Best.Oinf = @(X) OinfIfun(p(pOIidx),X(:,VarOinfI));
        fun_Best.DAinf = @(X) DAinfIfun(p(pDAIidx),X(:,VarDAinfI));
        fun_Best.Grect = @(X) Grectfun(p(pGidx),X(:,VarGrect));
    else
        if fval<fval_Best.I
            fval_Best.I = fval;
            idx_Best.I = [ifn];
            fun_Best.Oinf = @(X) OinfIfun(p(pOIidx),X(:,VarOinfI));
            fun_Best.DAinf = @(X) DAinfIfun(p(pDAIidx),X(:,VarDAinfI));
            fun_Best.Grect = @(X) Grectfun(p(pGidx),X(:,VarGrect));
        end
    end



    % plot figures
    if plot_flag
        tplot = 0:options.dt:options.tend;
        IChR2Onplot=iChR2On(p,xIp,tplot);
        [~,idxMplot]=max(abs(IChR2Onplot),[],2);
        idxM2plot=[[1:length(idxMplot)]',idxMplot];
        idxMplot=idxM2plot(:,1)+(idxM2plot(:,2)-1).*length(idxMplot);
        IpeakMplot=(IChR2Onplot(idxMplot));
        IssMplot=(IChR2Onplot(:,end));

        figure(fignr)
        if figpos_flag
            set(gcf,'position',figposition)
        end
        subplot(2,2,1)
        for ip=1:length(Il_valsp)
            idxx = find(xIp(:,2)==Il_valsp(ip));
            [xvals,idxvals] = sort(xIp(idxx,1));
            if length(xvals)>1
                plot(xvals,IssMplot(idxx(idxvals)),type(ifn),'color',Colors(ip),'DisplayName',sprintf('fit %5.2e',Il_valsp(ip)))
                hold on
                plot(xIss(idxx,1),yIss(idxx),'*','color',Colors(ip),'DisplayName',sprintf('data %5.2e',Il_valsp(ip)))
            end
        end
        hold off
        xlabel('Vm (mV)');
        ylabel('Iss')
        legend('show')

        subplot(2,2,2)
        for ip=1:length(Vm_valsp)
            idxx = find(xIp(:,1)==Vm_valsp(ip));
            [xvals,idxvals] = sort(xIp(idxx,2));
            if length(xvals)>1
                plot(xvals,IssMplot(idxx(idxvals)),type(ifn),'color',Colors(ip),'DisplayName',sprintf('fit %5.2e',Vm_valsp(ip)))
                hold on
                plot(xIss(idxx,2),yIss(idxx),'*','color',Colors(ip),'DisplayName',sprintf('data %5.2e',Vm_valsp(ip)))
            end
        end
        hold off
        xlabel('Il (W/m�)');
        legend('show')

        subplot(2,2,3)
        for ip=1:length(Il_valsp)
            idxx = find(xIp(:,2)==Il_valsp(ip));
            [xvals,idxvals] = sort(xIp(idxx,1));
            if length(xvals)>1
                plot(xvals,IpeakMplot(idxx(idxvals)),type(ifn),'color',Colors(ip),'DisplayName',sprintf('fit %5.2e',Il_valsp(ip)))
                hold on
                plot(xIp(idxx,1),yIp(idxx),'*','color',Colors(ip),'DisplayName',sprintf('data %5.2e',Il_valsp(ip)))
            end
        end
        hold off
        xlabel('Vm (mV)')
        ylabel('Ipeak')
        legend('show')

        subplot(2,2,4)
        for ip=1:length(Vm_valsp)
            idxx = find(xIp(:,1)==Vm_valsp(ip));
            [xvals,idxvals] = sort(xIp(idxx,2));
            if length(xvals)>1
                plot(xvals,IpeakMplot(idxx(idxvals)),type(ifn),'color',Colors(ip), 'DisplayName',sprintf('fit %5.2e',Vm_valsp(ip)))
                hold on
                plot(xIp(idxx,2),yIp(idxx),'*','color',Colors(ip),'DisplayName',sprintf('fit %5.2e',Vm_valsp(ip)))
            end
        end
        hold off
        xlabel('Il (W/m�)')
        legend('show')

        figure(fignr+1)
        if figpos_flag
            set(gcf,'position',figposition)
        end
        subplot(3,2,1)
        for ip=1:length(Il_valsp)

            idxx = find(xIp(:,2)==Il_valsp(ip));
            [xvals, idxvals] = sort(xIp(idxx,1));
            if length(xvals)>1
                plot(xvals,...
                    OinfIfun(p(pOIidx),xIp(idxx(idxvals),VarOinfI)),type(ifn),'color',Colors(ip),'DisplayName',sprintf('%5.2e',Il_valsp(ip)))
                hold on
            end
        end
        hold off
        xlabel('Il (W/m�)')
        ylabel('Oinf')
        legend('show')

        subplot(3,2,2)
        for ip=1:length(Vm_valsp)

            idxx = find(xIp(:,1)==Vm_valsp(ip),2);
            [xvals, idxvals] = sort(xIp(idxx,2));
            if length(xvals)>1
                plot(xvals,...
                    OinfIfun(p(pOIidx),xIp(idxx(idxvals),VarOinfI)),type(ifn),'color',Colors(ip),'DisplayName',sprintf('%5.2e',Vm_valsp(ip)))
                hold on
            end
        end
        hold off
        xlabel('Vm (mV)')
        legend('show')

        subplot(3,2,3)
        for ip=1:length(Il_valsp)
            idxx = find(xIp(:,2)==Il_valsp(ip));
            [xvals, idxvals] = sort(xIp(idxx,1));
            if length(xvals)>1
                plot(xvals,...
                    DAinfIfun(p(pDAIidx),xIp(idxx(idxvals),VarDAinfI)),type(ifn),'color',Colors(ip),'DisplayName',sprintf('%5.2e',Il_valsp(ip)))
                hold on
            end
        end
        hold off
        xlabel('Vm (mV)')
        ylabel('DAinf')
        legend('show')

        subplot(3,2,4)
        for ip=1:length(Vm_valsp)

            idxx = find(xIp(:,1)==Vm_valsp(ip));
            [xvals, idxvals] = sort(xIp(idxx,2));
            if length(xvals)>1
                plot(xvals,...
                    DAinfIfun(p(pDAIidx),xIp(idxx(idxvals),VarDAinfI)),type(ifn),'color',Colors(ip),'DisplayName',sprintf('%5.2e',Vm_valsp(ip)))
                hold on
            end
        end
        hold off
        xlabel('Vm (mV)')
        legend('show')

        subplot(3,2,5)
        for ip=1:length(Il_valsp)
            idxx = find(xIp(:,2)==Il_valsp(ip));
            [xvals, idxvals] = sort(xIp(idxx,1));
            if length(xvals)>1
                plot(xvals,...
                    Grectfun(p(pGidx),xIp(idxx(idxvals),VarGrect)),type(ifn),'color',Colors(ip),'DisplayName',sprintf('%5.2e',Il_valsp(ip)))
                hold on
            end
        end
        hold off
        xlabel('Vm (mV)')
        ylabel('Grect')
        legend('show')

        subplot(3,2,6)
        for ip=1:length(Vm_valsp)

            idxx = find(xIp(:,1)==Vm_valsp(ip));
            [xvals, idxvals] = sort(xIp(idxx,2));
            if length(xvals)>1
                plot(xvals,...
                    Grectfun(p(pGidx),xIp(idxx(idxvals),VarGrect)),type(ifn),'color',Colors(ip),'DisplayName',sprintf('%5.2e',Vm_valsp(ip)))
                hold on
            end
        end
        hold off
        xlabel('Vm (mV)')
        legend('show')
    end
end

    function OBJ=objfunIpeakIss(p,tX,tPeak,tSS,IChR2Onfun,options)
        to=0:options.dt:options.tend;
        IChR2On=IChR2Onfun(p,tX,to);
        [~,idxM]=max(abs(IChR2On),[],2);
        idxM2=[[1:length(idxM)]',idxM];
        idxM=idxM2(:,1)+(idxM2(:,2)-1).*length(idxM);
        IpeakM=(IChR2On(idxM));
        IssM=(IChR2On(:,end));
        OBJ_M=[options.wpeak.*(IpeakM-tPeak),options.wss.*(IssM-tSS)];
        OBJ=sqrt(mean(OBJ_M(:).^2));
    end
end