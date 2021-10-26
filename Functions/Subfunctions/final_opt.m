function [Out,fun_Best,idx_Best,fval_Best] = final_opt(selection,objMethod,solverMethod,OinfIfn,TauOIfn,...
    DAinfIfn,TauDAIfn,TauOVfn,TauDAVfn,Grectfn,combTauO,combTauDA,...
    param,fun,powera,powerb,Vm_valsp,Il_valsp,idx_Best,...
    xdata,ydata,Target_sp,Target_tr,Vm_valtr,options,display_flag,plot_flag,plotTauRecov_flag,TauRecov_flag,Colors,memorySave_flag)

if display_flag
    Collectionfunnames={'OinfI','TauOI','IinfI','TauII','TauOV','TauIV','Grect','combTauO','combTauDA'};
end

if any([length(OinfIfn),length(TauOIfn),length(DAinfIfn),length(TauDAIfn),length(TauOVfn),length(TauDAVfn),length(Grectfn)]-length(OinfIfn))
    error('length fun names do not match')
end


switch lower(selection)
    case 'all'
        for ifn=1:length(OinfIfn)

            [p,fval,OinfIfun,pOIidx,VarOinfI,pTOIidx,DAinfIfun,pDAIidx,VarDAinfI,pTDAIidx,pTOVidx,pTDAVidx,...
                Grectfun,pGidx,VarGrect,TauOVfun,TauOIfun,TauDAVfun,TauDAIfun,TauO,TauDA,X0,LB,UB,exitflag,output] =...
                runoptimization(fun,param,powera,powerb,Target_sp,Target_tr,objMethod,solverMethod,Vm_valsp,OinfIfn,TauOIfn,...
                DAinfIfn,TauDAIfn,TauOVfn,TauDAVfn,Grectfn,combTauO,combTauDA,ifn,ifn,ifn,options,...
                Collectionfunnames,TauRecov_flag,display_flag,plot_flag,Colors);

            Out.(selection).param(ifn).p = p;
            Out.(selection).param(ifn).pOI = p(pOIidx);
            Out.(selection).param(ifn).pOIidx = pOIidx;
            Out.(selection).param(ifn).pDAI = p(pDAIidx);
            Out.(selection).param(ifn).pDAIidx = pDAIidx;
            Out.(selection).param(ifn).pTOI = p(pTOIidx);
            Out.(selection).param(ifn).pTOIidx = pTOIidx;
            Out.(selection).param(ifn).pTDAI = p(pTDAIidx);
            Out.(selection).param(ifn).pTDAIidx = pTDAIidx;
            Out.(selection).param(ifn).pTOV = p(pTOVidx);
            Out.(selection).param(ifn).pTOVidx = pTOVidx;
            Out.(selection).param(ifn).pTDAV = p(pTDAVidx);
            Out.(selection).param(ifn).pTDAVidx = pTDAVidx;
            Out.(selection).param(ifn).pG = p(pGidx);
            Out.(selection).param(ifn).pGidx = pGidx;
            Out.(selection).param(ifn).fval = fval;
            Out.(selection).param(ifn).OIn = OinfIfn{ifn};
            Out.(selection).param(ifn).OIfun = func2str(OinfIfun);
            Out.(selection).param(ifn).DAIn = DAinfIfn{ifn};
            Out.(selection).param(ifn).DAIfun = func2str(DAinfIfun);
            Out.(selection).param(ifn).Gn = Grectfn{ifn};
            Out.(selection).param(ifn).Gfun = func2str(Grectfun);
            Out.(selection).param(ifn).TauOVn = TauOVfn{ifn};
            Out.(selection).param(ifn).TauOVfun = func2str(TauOVfun);
            Out.(selection).param(ifn).TauOIn = TauOIfn{ifn};
            Out.(selection).param(ifn).TauOIfun = func2str(TauOIfun);
            Out.(selection).param(ifn).TauDAVn = TauDAVfn{ifn};
            Out.(selection).param(ifn).TauDAVfun = func2str(TauDAVfun);
            Out.(selection).param(ifn).TauDAIn = TauDAIfn{ifn};
            Out.(selection).param(ifn).TauDAIfun = func2str(TauDAIfun);
            Out.(selection).param(ifn).combTauO = func2str(combTauO);
            Out.(selection).param(ifn).combTauDA = func2str(combTauDA);
            Out.(selection).objmethod = objMethod;
            Out.(selection).solvermethod = solverMethod;
            Out.(selection).powera = powera;
            Out.(selection).powerb = powerb;
            Out.(selection).output = output;
            Out.(selection).exitflag = exitflag;
            if ~memorySave_flag
                Out.(selection).param(ifn).X0 = X0;
                Out.(selection).param(ifn).LB = LB;
                Out.(selection).param(ifn).UB = UB;
                Out.(selection).options = options;
            end

            if ifn==1

                fval_Best.fopt = fval;
                idx_Best.fopt = [ifn];
                fun_Best.Oinf = @(X) OinfIfun(p(pOIidx),X(:,VarOinfI));
                fun_Best.DAinf = @(X) DAinfIfun(p(pDAIidx),X(:,VarDAinfI));
                fun_Best.Grect = @(X) Grectfun(p(pGidx),X(:,VarGrect));
                fun_Best.TauO = @(X) TauO(p,X);
                fun_Best.TauDA = @(X) TauDA(p,X);
            else
                if fval<fval_Best.fopt
                    fval_Best.fopt = fval;
                    idx_Best.fopt = [ifn];
                    fun_Best.Oinf = @(X) OinfIfun(p(pOIidx),X(:,VarOinfI));
                    fun_Best.DAinf = @(X) DAinfIfun(p(pDAIidx),X(:,VarDAinfI));
                    fun_Best.Grect = @(X) Grectfun(p(pGidx),X(:,VarGrect));
                    fun_Best.TauO = @(X) TauO(p,X);
                    fun_Best.TauDA = @(X) TauDA(p,X);
                end
            end
        end

    case 'best'
        ifI = idx_Best.I;
        ifTauO = idx_Best.TauO;
        ifTauDA = idx_Best.TauDA;

        [p,fval,OinfIfun,pOIidx,VarOinfI,pTOIidx,DAinfIfun,pDAIidx,VarDAinfI,pTDAIidx,pTOVidx,pTDAVidx,...
            Grectfun,pGidx,VarGrect,TauOVfun,TauOIfun,TauDAVfun,TauDAIfun,TauO,TauDA,X0,LB,UB,exitflag,output] =...
            runoptimization(fun,param,powera,powerb,Target_sp,Target_tr,objMethod,solverMethod,Vm_valsp,OinfIfn,TauOIfn,...
            DAinfIfn,TauDAIfn,TauOVfn,TauDAVfn,Grectfn,combTauO,combTauDA,ifI,ifTauO,ifTauDA,options,...
            Collectionfunnames,TauRecov_flag,display_flag,plot_flag,Colors);
        if isempty(p)
            warning('empty p final opt')
            p =X0;
            fval = inf;
        end

        Out.(selection).p = p;
        Out.(selection).pOI = p(pOIidx);
        Out.(selection).pOIidx = pOIidx;
        Out.(selection).pDAI = p(pDAIidx);
        Out.(selection).pDAIidx = pDAIidx;
        Out.(selection).pTOI = p(pTOIidx);
        Out.(selection).pTOIidx = pTOIidx;
        Out.(selection).pTDAI = p(pTDAIidx);
        Out.(selection).pTDAIidx = pTDAIidx;
        Out.(selection).pTOV = p(pTOVidx);
        Out.(selection).pTOVidx = pTOVidx;
        Out.(selection).pTDAV = p(pTDAVidx);
        Out.(selection).pTDAVidx = pTDAVidx;
        Out.(selection).pG = p(pGidx);
        Out.(selection).pGidx = pGidx;
        Out.(selection).fval = fval;
        Out.(selection).OIn = OinfIfn{ifI};
        Out.(selection).OIfun = func2str(OinfIfun);
        Out.(selection).DAIn = DAinfIfn{ifI};
        Out.(selection).DAIfun = func2str(DAinfIfun);
        Out.(selection).Gn = Grectfn{ifI};
        Out.(selection).Gfun = func2str(Grectfun);
        Out.(selection).TauOVn = TauOVfn{ifTauO};
        Out.(selection).TauOVfun = func2str(TauOVfun);
        Out.(selection).TauOIn = TauOIfn{ifTauO};
        Out.(selection).TauOIfun = func2str(TauOIfun);
        Out.(selection).TauDAVn = TauDAVfn{ifTauDA};
        Out.(selection).TauDAVfun = func2str(TauDAVfun);
        Out.(selection).TauDAIn = TauDAIfn{ifTauDA};
        Out.(selection).TauDAIfun = func2str(TauDAIfun);
        Out.(selection).combTauO = func2str(combTauO);
        Out.(selection).combTauDA = func2str(combTauDA);
        Out.(selection).objmethod = objMethod;
        Out.(selection).solvermethod = solverMethod;
        Out.(selection).powera = powera;
        Out.(selection).powerb = powerb;
        Out.(selection).output = output;
        Out.(selection).exitflag = exitflag;
        if ~memorySave_flag
            Out.(selection).X0 = X0;
            Out.(selection).LB = LB;
            Out.(selection).UB = UB;
            Out.(selection).options = options;
        end

        fval_Best.fopt = fval;
        idx_Best.fopt = [ifn];
        fun_Best.Oinf = @(X) OinfIfun(p(pOIidx),X(:,VarOinfI));
        fun_Best.DAinf = @(X) DAinfIfun(p(pDAIidx),X(:,VarDAinfI));
        fun_Best.Grect = @(X) Grectfun(p(pGidx),X(:,VarGrect));
        fun_Best.TauO = @(X) TauO(p,X);
        fun_Best.TauDA = @(X) TauDA(p,X);

    otherwise
        error ('false input')
end


    function [p,fval,OinfIfun,pOIidx,VarOinfI,pTOIidx,DAinfIfun,pDAIidx,VarDAinfI,pTDAIidx,pTOVidx,pTDAVidx,...
            Grectfun,pGidx,VarGrect,TauOVfun,TauOIfun,TauDAVfun,TauDAIfun,TauO,TauDA,X0,LB,UB,exitflag,output] =...
            runoptimization(fun,param,powera,powerb,Target_sp,Target_tr,objMethod,solverMethod,Vm_valsp,OinfIfn,TauOIfn,DAinfIfn,TauDAIfn,TauOVfn,TauDAVfn,Grectfn,combTauO,combTauDA,ifI,ifTauO,ifTauDA,options,...
            Collectionfunnames,TauRecov_flag,display_flag,plot_flag,Colors)

        incl_tr_flag = options.incl_tr_flag;
        if display_flag
            Selectedfunnames={OinfIfn{ifI},TauOIfn{ifTauO},DAinfIfn{ifI},TauDAIfn{ifTauDA},TauOVfn{ifTauO},TauDAVfn{ifTauDA},...
                Grectfn{ifI},func2str(combTauO),func2str(combTauDA)};
            Stringtoprint={Collectionfunnames{:};Selectedfunnames{:}};
            fprintf('\n selected functions: %s = %s \n',Stringtoprint{:})
        end

        OinfIfun = fun.(OinfIfn{ifI}).fun;
        Pidx(1) = fun.(OinfIfn{ifI}).nrpara;
        pOIidx = 1:Pidx(1);
        X0OinfI = param.I.param(ifI).pOI;
        LBOinfI_init = fun.(OinfIfn{ifI}).LB;
        UBOinfI_init = fun.(OinfIfn{ifI}).UB;
        VarOinfI = fun.(OinfIfn{ifI}).var;

        TauOIfun = fun.(TauOIfn{ifTauO}).fun;
        Pidx(2) = fun.(TauOIfn{ifTauO}).nrpara+Pidx(1);
        pTOIidx = Pidx(1)+1:Pidx(2);
        X0TauOI = param.TauO.param(ifTauO).pI;
        LBTauOI_init = fun.(TauOIfn{ifTauO}).LB;
        UBTauOI_init = fun.(TauOIfn{ifTauO}).UB;
        VarTauOI = fun.(TauOIfn{ifTauO}).var;

        DAinfIfun = fun.(DAinfIfn{ifI}).fun;
        Pidx(3) = fun.(DAinfIfn{ifI}).nrpara+Pidx(2);
        pDAIidx = Pidx(2)+1:Pidx(3);
        X0DAinfI = param.I.param(ifI).pDAI;
        LBDAinfI_init = fun.(DAinfIfn{ifI}).LB;
        UBDAinfI_init = fun.(DAinfIfn{ifI}).UB;
        VarDAinfI = fun.(DAinfIfn{ifI}).var;

        TauDAIfun=fun.(TauDAIfn{ifTauDA}).fun;
        Pidx(4)=fun.(TauDAIfn{ifTauDA}).nrpara+Pidx(3);
        pTDAIidx = Pidx(3)+1:Pidx(4);
        X0TauDAI=param.TauDA.param(ifTauDA).pI;
        LBTauDAI_init = fun.(TauDAIfn{ifTauDA}).LB;
        UBTauDAI_init = fun.(TauDAIfn{ifTauDA}).UB;
        VarTauDAI=fun.(TauDAIfn{ifTauDA}).var;

        TauOVfun=fun.(TauOVfn{ifTauO}).fun;
        Pidx(5)=fun.(TauOVfn{ifTauO}).nrpara+Pidx(4);
        pTOVidx = Pidx(4)+1:Pidx(5);
        X0TauOV=param.TauO.param(ifTauO).pV;
        LBTauOV_init = fun.(TauOVfn{ifTauO}).LB;
        UBTauOV_init = fun.(TauOVfn{ifTauO}).UB;
        VarTauOV=fun.(TauOVfn{ifTauO}).var;

        TauDAVfun=fun.(TauDAVfn{ifTauDA}).fun;
        Pidx(6)=fun.(TauDAVfn{ifTauDA}).nrpara+Pidx(5);
        pTDAVidx = Pidx(5)+1:Pidx(6);
        X0TauDAV=param.TauDA.param(ifTauDA).pV;
        LBTauDAV_init = fun.(TauDAVfn{ifTauDA}).LB;
        UBTauDAV_init = fun.(TauDAVfn{ifTauDA}).UB;
        VarTauDAV=fun.(TauDAVfn{ifTauDA}).var;

        Grectfun=fun.(Grectfn{ifI}).fun;
        Pidx(7)=Pidx(6)+fun.(Grectfn{ifI}).nrpara;
        pGidx = Pidx(6)+1:Pidx(7);
        X0Grect=param.I.param(ifI).pG;
        LBGrect_init = fun.(Grectfn{ifI}).LB;
        UBGrect_init = fun.(Grectfn{ifI}).UB;
        VarGrect=fun.(Grectfn{ifI}).var;

        TauO = @(p,X) combTauO(TauOVfun(p(pTOVidx),X(:,VarTauOV)),...
            TauOIfun(p(pTOIidx),X(:,VarTauOI)));
        TauDA = @(p,X) combTauDA(TauDAVfun(p(pTDAVidx),X(:,VarTauDAV)),...
            TauDAIfun(p(pTDAIidx),X(:,VarTauDAI)));

        yfunx = @(t,taux,xinf,x0) xinf+(x0-xinf).*exp(-(t.*double(t>=0))./taux);

        Oton = @(p,X,t) yfunx(t,TauO(p,X),OinfIfun(p(pOIidx),X(:,VarOinfI)),0);
        DAton = @(p,X,t) yfunx(t,TauDA(p,X),DAinfIfun(p(pDAIidx),X(:,VarDAinfI)),1);

        Otoff = @(p,X,t,dton) yfunx(t,TauO(p,[X(:,1),zeros(size(X(:,1)))]),0,Oton(p,X,dton));
        DAtoff = @(p,X,t,dton) yfunx(t,TauDA(p,[X(:,1),zeros(size(X(:,1)))]),1,DAton(p,X,dton));

        Oton_2p = @(p,X,t,dtoff,dton) yfunx(t,TauO(p,X),OinfIfun(p(pOIidx),X(:,VarOinfI)),Otoff(p,X,dtoff,dton));
        DAton_2p = @(p,X,t,dtoff,dton) yfunx(t,TauDA(p,X),DAinfIfun(p(pDAIidx),X(:,VarDAinfI)),DAtoff(p,X,dtoff,dton));

        Otoff_2p = @(p,X,t,dton2,dtoff1,dton1) yfunx(t,TauO(p,[X(:,1),zeros(size(X(:,1)))]),0,Oton_2p(p,X,dton2,dtoff1,dton1));
        DAtoff_2p = @(p,X,t,dton2,dtoff1,dton1) yfunx(t,TauDA(p,[X(:,1),zeros(size(X(:,1)))]),1,DAton_2p(p,X,dton2,dtoff1,dton1));

        iChR2On=@(p,X,t) Grectfun(p(pGidx),X(:,VarGrect)).*(Oton(p,X,t)).^powera.*(DAton(p,X,t)).^powerb;
        iChR2Off=@(p,X,t,dton) Grectfun(p(pGidx),X(:,VarGrect)).*Otoff(p,X,t,dton).^powera.*(DAtoff(p,X,t,dton)).^powerb;

        iChR2On_2p=@(p,X,t,dtoff,dton) Grectfun(p(pGidx),X(:,VarGrect)).*(Oton_2p(p,X,t,dtoff,dton)).^powera.*(DAton_2p(p,X,t,dtoff,dton)).^powerb;
        iChR2Off_2p=@(p,X,t,dton2,dtoff1,dton1) Grectfun(p(pGidx),X(:,VarGrect)).*Otoff_2p(p,X,t,dton2,dtoff1,dton1).^powera.*(DAtoff_2p(p,X,t,dton2,dtoff1,dton1)).^powerb;

        iChR2fun_sp = @(p,X,t,ton1,toff1) iChR2On(p,X,t-ton1).*double(t>=ton1&t<toff1)+iChR2Off(p,X,t-toff1,toff1-ton1).*double(t>=toff1);

        iChR2fun_2p = @(p,X,t,ton1,toff1,ton2,toff2) iChR2On(p,X,max(t-ton1,0)).*double(t>=ton1&t<toff1)+iChR2Off(p,X,max(t-toff1,0),toff1-ton1).*double(t>=toff1&t<ton2) +...
            iChR2On_2p(p,X,max(t-ton2,0),ton2-toff1,toff1-ton1).*(t>=ton2&t<toff2) + iChR2Off_2p(p,X,max(t-toff2,0),toff2-ton2,ton2-toff1,toff1-ton1).*double(t>=toff2);

        % range based on previous fits
        LB_init = [LBOinfI_init,LBTauOI_init,LBDAinfI_init,LBTauDAI_init,LBTauOV_init,LBTauDAV_init,LBGrect_init];
        UB_init = [UBOinfI_init,UBTauOI_init,UBDAinfI_init,UBTauDAI_init,UBTauOV_init,UBTauDAV_init,UBGrect_init];

        X0 = [X0OinfI,X0TauOI,X0DAinfI,X0TauDAI,X0TauOV,X0TauDAV,X0Grect];
        switch lower(options.ftype)
            case 'scale'
                LB = X0.*options.fLB;
                UB = X0.*options.fUB;
            case 'add'
                LB = X0-options.fLB;
                UB = X0+options.fUB;
            case 'frac_initb'
                rLBUB = UB_init-LB_init;
                LB = X0-options.fLB.*rLBUB;
                UB = X0+options.fUB.*rLBUB;
            otherwise
                error('wrong scale method LB , UB')
        end

        UB2=UB;
        UB(UB2<LB)=LB(UB2<LB);
        LB(UB2<LB)=UB2(UB2<LB);

        UB = min(UB_init,UB);
        LB = max(LB_init,LB);

        tV_sp = [Target_sp(:).Vm];  tIl_sp = [Target_sp(:).Il];
        ton1_sp = [Target_sp(:).OSpstart];  toff1_sp = ton1_sp + [Target_sp(:).OSpd];
        if TauRecov_flag
            ton1_tr = [Target_tr(:).Pstart];    toff1_tr = [Target_tr(:).PD]+ton1_tr;
            tVr = [Target_tr(:).Vm]; tIr = [Target_tr(:).Il];
        else

            ton1_tr = [];    toff1_tr = [];
            tVr = []; tIr = [];
        end

        switch lower(objMethod)
            case 'currenttraces'

                OBJ = @(p) objfun_CT(p,tV_sp,tI_sp,ton1_sp,toff1_sp,Target_sp,...
                    iChR2fun_sp,ton1_tr,toff1_tr,tVr,tIr,iChR2fun_2p,Target_tr,incl_tr_flag,options.CT);

            case 'features'

                tTauOn = [Target_sp(:).TauOn];  tTauOff = [Target_sp(:).TauOff];
                tTauInact = [Target_sp(:).TauInact];    tTauRecov = [Target_tr(:).TauRecov];
                tIpeak = [Target_sp(:).Ipeak];  tIss = [Target_sp(:).Iss];
                tIratio = [Target_sp(:).Iratio];

                pos_nantargets=[isnan(tTauOn'),isnan(tTauInact'),isnan(tTauOff'),isnan(tIpeak'),isnan(tIss'),isnan(tIratio')];

                if isfield(Target_sp,'time')
                    tend_sp = cellfun(@(t) t(end),{Target_sp(:).time})';
                else
                    tend_sp = options.FB.tend_sp;
                    if length(tend_sp)~=length(tV_sp)
                        if length(tend_sp) ==1
                            tend_sp = ones(length(tV_sp),1)*tend_sp;
                        else
                            error('false length(tend_sp)')
                        end
                    end
                end

                if TauRecov_flag
                    if isfield(Target_tr,'tend')
                        tend_tr = vertcat(Target_tr(:).tend);
                    else
                        tend_tr = options.FB.tend_tr;
                        if length(tend_tr)~=length(tV_tr)
                            if length(tend_tr) ==1
                                tend_tr= ones(length(tV_tr),length(options.FB.intervals))*tend_tr;
                            else
                                error('false length(tend_sp)')
                            end
                        end
                    end
                end

                CF = 1-log(1/(1-nanmedian(tIratio)));
                if CF<0.1
                    CF = 0.1;
                    warning('CF for Tau Recov fixed to 0.1')
                end
                feat_extr_settings = {'ssReached_flag',options.FB.ssReached_flag,'segTauInact_startpoint',options.FB.segTauInact_startpoint};
                if isfield(options.FB,'feat_extr_settings')
                    feat_extr_settings = horzcat(feat_extr_settings,options.FB.feat_extr_settings);
                end

                OBJ = @(p) objfun_FB(p,tV_sp,tIl_sp,ton1_sp,toff1_sp,ton1_tr,toff1_tr,tVr,tIr,tTauOn,tTauOff,tTauInact,tTauRecov,...
                    tIpeak,tIss,tIratio,iChR2fun_sp,iChR2fun_2p,TauDA,powerb,tend_sp,tend_tr,incl_tr_flag,CF,pos_nantargets,feat_extr_settings,options.FB);

            otherwise
                error('false method')
        end

        %preparing evaluation values for nonlincon
        XO = options.I_nonlcon;
        XDA = options.I_nonlcon;
        XDAend = tIl_sp; if isrow(XDAend); XDAend = XDAend'; end
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
            XDAend = [tV_sp(:),tIl_sp(:)];
            XDA = XDA(:)';
            Vvals = [XV1(:);XV2(:)];
            [Vvals,XDA] = meshgrid(Vvals,XDA);
            XDA = [Vvals(:),XDA(:)];
        elseif numel(VarDAinfI)==1 &&  size(XDA,2)==2
            error('check input nonlincon both IV but funtion single column required')
        elseif numel(VarDAinfI)==1 &&  isrow(XDA)
            XDA = XDA(:);
        end

        nonlcon = @(p) nonlincon(p,powera,powerb,OinfIfun,pOIidx,DAinfIfun,pDAIidx,Grectfun,pGidx,XO,XO0,XDA,XDA0,XV1,XV2,...
            options.extrac_nonlcon,options.extraceq_nonlcon,TauO(p,[tV_sp',zeros(size(tV_sp'))]),TauDA(p,[tV_sp',zeros(size(tV_sp'))]),XDAend);
        outFun = @(x,optimValues,state) outfun(x,optimValues,state,...
            plot_flag,OinfIfun,pOIidx,VarOinfI,DAinfIfun,pDAIidx,VarDAinfI,TauO,TauDA,Grectfun,...
            pGidx,VarGrect,Vm_valsp,Target_sp,options.Outfun,powera,powerb,Colors,0,tVr,Target_tr);
        switch lower(solverMethod)
            case 'ms'
                disp('start ms finalopt')
                opt = optimoptions(options.ms.optimoptions{:},'OutputFcn',outFun);
                problem=createOptimProblem('fmincon',...
                    'objective',OBJ,...
                    'nonlcon',nonlcon,...
                    'x0',X0,...
                    'lb',LB,...
                    'ub',UB,...
                    'options',opt);
                ms = MultiStart(options.ms.msopt{:});
                [p,fval,exitflag,output]= run(ms,problem,options.ms.msnr);

            case 'psobc'
                disp('start psoBC finalopt')
                problem = struct();
                problem = options.psoBC; % add simulation options
                problem.options.InitialPopulation = [X0];
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

                [p,fval,exitflag,output] = pso(problem);
            case 'ga'
                disp('start ga finalopt')
                problem = struct();
                problem = options.ga; % add simulation options
                problem.options.InitialPopulation = [X0];
                problem.options.creationFunction = 'gacreationnonlinearfeasible';
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

                [p,fval,exitflag,output] = ga(problem);
            otherwise
                error('false input')
        end
    end
end