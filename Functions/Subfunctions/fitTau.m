function [Out,fun_Best,idx_Best,fval_Best] = fitTau(xdata,ydata,Ttype,fun,...
    TauIfun_names,TauVfun_names,combTau,options,plot_flag,IlTau_val,Vm_valsp,Colors,type,memorySave_flag,plotTau_flag,varargin)

nonlcon_flag = 0;

TauDAflag = double(strcmpi(Ttype,'TauDA'));
fignr = get(gcf,'number')+(1-TauDAflag)*(1-plotTau_flag);

% check if TauIfn_names and TauVfun_names same length
if length(TauIfun_names) ~= length(TauVfun_names)
    error('TauIfun_names has to have the same length as TauVfun_names')
end

if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'nonlcon'))
            nonlcon_flag = varargin{find(strcmpi(varargin,'nonlcon'))+1};
            if nonlcon_flag
                TauDA = varargin{find(strcmpi(varargin,'TauDA'))+1};
                DAinf = varargin{find(strcmpi(varargin,'DAinf'))+1};
                powera = varargin{find(strcmpi(varargin,'powera'))+1};
                powerb = varargin{find(strcmpi(varargin,'powerb'))+1};
            end
        end
    end
else
    error('varargin wrong size')
end

for ifn = 1:length(TauIfun_names)

    TauIfun=fun.(TauIfun_names{ifn}).fun;
    Pidx(1)=fun.(TauIfun_names{ifn}).nrpara;
    pTIidx = 1:Pidx(1);
    LBTauI=fun.(TauIfun_names{ifn}).LB;
    UBTauI=fun.(TauIfun_names{ifn}).UB;
    X0TauI=fun.(TauIfun_names{ifn}).X0;
    VarTauI=fun.(TauIfun_names{ifn}).var;

    TauVfun=fun.(TauVfun_names{ifn}).fun;
    Pidx(2)=fun.(TauVfun_names{ifn}).nrpara+Pidx(1);
    pTVidx = Pidx(1)+1:Pidx(2);
    LBTauV=fun.(TauVfun_names{ifn}).LB;
    UBTauV=fun.(TauVfun_names{ifn}).UB;
    X0TauV=fun.(TauVfun_names{ifn}).X0;
    VarTauV=fun.(TauVfun_names{ifn}).var;

    X0=[X0TauI,X0TauV];
    LB=[LBTauI,LBTauV];
    UB=[UBTauI,UBTauV];

    % combination function requirec V dependence first
    Taufun = @(p,X) combTau(TauVfun(p(pTVidx),X(:,VarTauV)),...
        TauIfun(p(pTIidx),X(:,VarTauI)));


    switch options.solverMethod
        case 'lsqcurvefit'

            opt=optimoptions('lsqcurvefit',options.optimoptions{:});

            if nonlcon_flag && ~TauDAflag
                xdata_nlc = [xdata(:,1),zeros(size(xdata(:,2)))];
                DAinf_nlc = DAinf(xdata);
                TauDAoff = TauDA(xdata_nlc);
                nonlcon = @(p) nonlinconTau(p,powera,powerb,Taufun,xdata_nlc,DAinf_nlc,TauDAoff);
                problem=createOptimProblem('lsqcurvefit',...
                    'objective',Taufun,...
                    'xdata',xdata,'ydata',ydata,...
                    'nonlcon',nonlcon,...
                    'x0',X0,...
                    'lb',LB,...
                    'ub',UB, 'options',opt);
            else
                problem=createOptimProblem('lsqcurvefit',...
                    'objective',Taufun,...
                    'xdata',xdata,'ydata',ydata,...
                    'x0',X0,...
                    'lb',LB,...
                    'ub',UB, 'options',opt);
            end
        case 'fmincon'
            opt=optimoptions('fmincon','Algorithm','sqp',options.optimoptions{:});

            idxIlz = xdata(:,2)==0;
            xdata_Ilnz = xdata(~idxIlz,:);
            ydata_Ilnz = ydata(~idxIlz,:);
            xdata_Ilz = xdata(idxIlz,:);
            ydata_Ilz = ydata(idxIlz,:);
            wIlnz = 1;
            wIlz = sum(~idxIlz)./sum(idxIlz);

            OBJ = @(p)  objfun_Tau(p,xdata_Ilnz,ydata_Ilnz,xdata_Ilz,ydata_Ilz,wIlnz,wIlz,Taufun);
            if nonlcon_flag && ~TauDAflag
                xdata_nlc = [xdata(:,1),zeros(size(xdata(:,2)))];
                DAinf_nlc = DAinf(xdata);
                TauDAoff = TauDA(xdata_nlc);
                nonlcon = @(p) nonlinconTau(p,powera,powerb,Taufun,xdata_nlc,DAinf_nlc,TauDAoff);
                problem=createOptimProblem('fmincon',...
                    'objective',OBJ,...
                    'nonlcon',nonlcon,...
                    'x0',X0,...
                    'lb',LB,...
                    'ub',UB, 'options',opt);
            else
                problem=createOptimProblem('fmincon',...
                    'objective',OBJ,...
                    'x0',X0,...
                    'lb',LB,...
                    'ub',UB, 'options',opt);
            end
        otherwise
            error('not right method')
    end

    ms = MultiStart(options.msopt{:});
    [p,fval,exitflag,output,solutions]= run(ms,problem,options.msnr);
    Error=fval/mean(ydata).^2;
    R2=1-fval./sum((ydata-mean(ydata)).^2);

    % store in output

    Out.(Ttype).param(ifn).p = p;
    Out.(Ttype).param(ifn).pI= p(pTIidx);
    Out.(Ttype).param(ifn).pIidx = pTIidx;
    Out.(Ttype).param(ifn).pV= p(pTVidx);
    Out.(Ttype).param(ifn).pVidx = pTVidx;
    Out.(Ttype).param(ifn).fval = fval;
    Out.(Ttype).param(ifn).error = Error;
    Out.(Ttype).param(ifn).R2 = R2;
    Out.(Ttype).param(ifn).TauVn = TauVfun_names{ifn};
    Out.(Ttype).param(ifn).TauVfun = func2str(TauVfun);
    Out.(Ttype).param(ifn).TauIn = TauIfun_names{ifn};
    Out.(Ttype).param(ifn).TauIfun = func2str(TauIfun);
    Out.(Ttype).param(ifn).combTau = func2str(combTau);
    if ~memorySave_flag
        Out.(Ttype).param(ifn).X0 = X0;
        Out.(Ttype).param(ifn).LB = LB;
        Out.(Ttype).param(ifn).UB = UB;
        Out.(Ttype).xdata = xdata;
        Out.(Ttype).ydata = ydata;
        Out.(Ttype).options = options;
    end

    if ifn==1
        fval_Best.(Ttype) = fval;
        idx_Best.(Ttype) = ifn;
        fun_Best.(Ttype) = @(X) Taufun(p,X);
    else
        if fval<fval_Best.(Ttype)
            fval_Best.(Ttype) = fval;
            idx_Best.(Ttype) = [ifn];
            fun_Best.(Ttype) = @(X) Taufun(p,X);
        end
    end


    if plot_flag
        cidx = 0;
        figure(fignr)
        subplot(2,2,1+2*TauDAflag)
        for ip=1:length(IlTau_val)
            idxx = find(xdata(:,2)==IlTau_val(ip));
            [xvals,idxvals] = sort(xdata(idxx,1));
            if length(xvals)>1
                cidx=cidx+1;
                hold on
                plot(xvals,Taufun(p,xdata(idxx(idxvals),:)),type(ifn),'color',Colors(cidx),'DisplayName',sprintf('fit %5.2e',IlTau_val(ip)))

            end
        end
        hold off

        cidx = 0;
        subplot(2,2,2+2*TauDAflag)
        for ip=1:length(Vm_valsp)
            idxx = find(xdata(:,1)==Vm_valsp(ip));
            [xvals,idxvals] = sort(xdata(idxx,2));
            if length(xvals)>1
                cidx=cidx+1;
                hold on
                plot(xvals,Taufun(p,xdata(idxx(idxvals),:)),type(ifn),'color',Colors(cidx),'DisplayName',sprintf('fit %5.2e',Vm_valsp(ip)))
            end
        end
        hold off
    end
end
end