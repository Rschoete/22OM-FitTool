function OBJ=objfun_FB(p,tV,tI,ton1_sp,toff1_sp,ton1_tr,toff1_tr,tVr,tIr,tTauOn,tTauOff,tTauInact,tTauRecov,...
    tIPeak,tIss,tRatio,iChR2fun_sp,iChR2fun_2p,TauDAfun,powerb,tend_sp,tend_tr,incl_tr_flag,CF,pos_nantargets,feat_extr_settings,options)
try
    if incl_tr_flag
        % Taurecov
        switch lower(options.method_tr)
            case 'estimpowerb'
                if powerb~=1
                TauRecov=estimTauRecov(TauDAfun(p,[tVr',zeros(size(tVr))']),tVr,powerb,options.trestim);
                else
                    TauRecov = TauDAfun(p,[tVr',zeros(size(tVr))']);
                end

            case 'noestimpowerb'
                TauRecov = TauDAfun(p,[tVr',zeros(size(tVr))']);
            case 'currenttraces'

                X_tr = [tVr',tIr'];
                intervals = options.intervals;
                TauRecov  = zeros(length(tVr),1);
                for itVr = 1:length(tVr)
                    Peak1 = zeros(length(intervals),1);
                    Peak2 = zeros(length(intervals),1);
                    Ratio = zeros(length(intervals),1);
                    for iInt = 1:length(intervals)
                        to_tr = 0:options.dt_tr:tend_tr(itVr,iInt);
                        ton2_tr = toff1_tr(itVr)+intervals(iInt);
                        toff2_tr = ton2_tr+toff1_tr(itVr)-ton1_tr(itVr);
                        iChR2 = iChR2fun_2p(p,X_tr(itVr,:),to_tr,ton1_tr(itVr),toff1_tr(itVr),ton2_tr,toff2_tr);
                        loc1=to_tr>=ton1_tr(itVr) & to_tr<=(toff1_tr(itVr)+1/2*intervals(iInt));
                        loc2=to_tr>=ton2_tr & to_tr<=toff2_tr+1/2*intervals(iInt);
                        Peak1(iInt)=max(abs(iChR2(loc1)));
                        %Peak2(i)=max(abs(Out(i).iChR2(loc2)));
                        try
                            Peaks2=findpeaks(abs(iChR2(loc2)));
                            Peak2(iInt)=Peaks2(1);
                        catch
                            Peak2(iInt) = 0;
                        end

                        Ratio(iInt)=Peak2(iInt)/Peak1(iInt);

                    end

                    xx = [0:options.dxx:intervals(end)];
                    RATIO=spline(intervals,Ratio,xx);
                    [~,idx]=min(abs(RATIO-(1-exp(-1))));
                    if any(Ratio<=(1-exp(-1))) && any(Ratio>=(1-exp(-1)))
                        TauRecov(itVr,1)=xx(idx(1));
                    else
                        if powerb~=1
                            TauRecov(itVr,1)=estimTauRecov(TauDAfun(p,[tVr',zeros(size(tVr))']),tVr,powerb,options.trestim)/CF;
                        else
                            TauRecov(itVr,1) = TauDAfun(p,[tVr(itVr),0])/CF;
                        end
                    end
                end

            case 'currenttraces_otherdef'

                X_tr = [tVr',tIr'];
                intervals = options.intervals;
                TauRecov  = zeros(length(tVr),1);
                for itVr = 1:length(tVr)
                    Peak1 = zeros(length(intervals),1);
                    Peak2 = zeros(length(intervals),1);
                    Ioff1 = zeros(length(intervals),1);
                    Ratio = zeros(length(intervals),1);
                    for iInt = 1:length(intervals)
                        to_tr = 0:options.dt_tr:tend_tr(itVr,iInt);
                        ton2_tr = toff1_tr(itVr)+intervals(iInt);
                        toff2_tr = ton2_tr+toff1_tr(itVr)-ton1_tr(itVr);
                        iChR2 = iChR2fun_2p(p,X_tr(itVr,:),to_tr,ton1_tr(itVr),toff1_tr(itVr),ton2_tr,toff2_tr);
                        loc1 = to_tr>=ton1_tr(itVr) & to_tr<=(toff1_tr(itVr)+1/2*intervals(iInt));
                        loc2=to_tr>=ton2_tr & to_tr<=toff2_tr+1/2*intervals(iInt);
                        Peak1(iInt)=max(abs(iChR2(loc1)));
                        loc_Ioff1 = to_tr <=toff1_tr(itVr);
                        Ioff1_intm = iChR2(loc_Ioff1);
                        Ioff1(iInt)=abs(Ioff1_intm(end)); Ioff1_intm = [];
                        %Peak2(i)=max(abs(Out(i).iChR2(loc2)));
                        try
                            Peaks2=findpeaks(abs(iChR2(loc2)));
                            Peak2(iInt)=Peaks2(1);
                        catch
                            Peak2(iInt) = 0;
                        end

                        Ratio(iInt)=(Peak2(iInt)-Ioff1(iInt))/(Peak1(iInt)-Ioff1(iInt));
                    end

                    xx = [0:options.dxx:intervals(end)];
                    RATIO=spline(intervals,Ratio,xx);
                    [~,idx]=min(abs(RATIO-(1-exp(-1))));
                    if any(Ratio<=(1-exp(-1))) && any(Ratio>=(1-exp(-1)))
                        TauRecov(itVr,1)=xx(idx(1));
                    else
                        %TauRecov(itVr,1)=nan;
                        if powerb~=1
                            TauRecov(itVr,1)=estimTauRecov(TauDAfun(p,[tVr',zeros(size(tVr))']),tVr,powerb,options.trestim);
                        else
                            TauRecov(itVr,1) = TauDAfun(p,[tVr(itVr),0]);
                        end
                    end

                end

            otherwise
                error('false method_tr')
        end
    end

    %single pulse
    X=[tV',tI']; %concatenate target holding voltage and pulse intensity 
    RMSMat=zeros(length(tV),6); % init cost matrix 
    for io=1:length(tV)
        to=0:options.dt:tend_sp(io); % define time points for simulation io
        iChR2_sp = iChR2fun_sp(p,X(io,:),to,ton1_sp(io),toff1_sp(io)); % calculate optocurrent under conditions X(io,:) for a voltage clamp at V=tV(io) and a pulse intensity of I=tI(io)
        Features=Extract_feat(to,iChR2_sp,toff1_sp(io)-ton1_sp(io),ton1_sp(io),options.nrr_featextr,options.plot_featextr,feat_extr_settings); %Extract features of calculated opto current
        
        %calculate weighted RMS error and store in cost matrix
        RMSMat(io,:)=options.w(io,:).*[(tTauOn(io)-Features.TauOn),(tTauInact(io)-Features.TauInact),(tTauOff(io)-Features.TauOff),...
            (tIPeak(io)-Features.Ipeak),(tIss(io)-Features.Iss),(tRatio(io)-Features.Iratio)];
    end
    %RMSMat=[RMSMat(:);(tTauRecov'-TauRecov)./tTauRecov'];
    RMSMat(pos_nantargets) = [];
    if incl_tr_flag
        RMSMat=[RMSMat(:);options.wtr.*(tTauRecov'-TauRecov)];
    end

    RMSMat(isnan(RMSMat))=10000;
    OBJ=sqrt(mean((RMSMat).^2));

catch e
    fprintf(1,'\nThe identifier was:\n%s',e.identifier);
    fprintf(1,'There was an error! The message was:\n%s',e.message)
    OBJ=10000;
    if any(size(options.w)-[length(tV),6])
        error('options.w false size')
    end
end
end