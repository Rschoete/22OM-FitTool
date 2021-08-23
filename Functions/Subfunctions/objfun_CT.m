function OBJ = objfun_CT(p,tV,tI,ton1_sp,toff1_sp,Target_sp,iChR2fun_sp,ton1_tr,toff1_tr,tVr,tIr,iChR2fun_2p,Target_tr,incl_tr_flag,options)
try
    %single pulse
    X=[tV',tI'];
    X_tr = [tVr',tIr'];
    intervals = options.intervals;
    switch options.method_CT
        case 'normalized'
            diffsq_mat = 0;
            nsamples = 0;
            iChR2_in = 0;
            for io=1:length(tV)

                iChR2_sp = iChR2fun_sp(p,X(io,:),Target_sp(io).time,ton1_sp(io),toff1_sp(io));
                iChR2_peakval = iChR2_sp(max(abs(iChR2_sp))==abs(iChR2_sp));

                iChR2_measured = Target_sp(io).I;
                iChR2_m_peakval = Target_sp(io).Ipeak;

                if iChR2_m_peakval~=0 && iChR2_peakval~=0
                    diffsq_mat = diffsq_mat+sum((iChR2_sp/iChR2_peakval-iChR2_measured/iChR2_m_peakval).^2);
                    nsamples = nsamples+length(iChR2_sp);
                    iChR2_in = iChR2_in + sum((iChR2_measured/iChR2_m_peakval).^2);
                end
            end
            if incl_tr_flag
                for itVr=1:length(X_tr)
                    for iInt=1:length(intervals)
                        ton2_tr = toff1_tr(itVR)+intervals(iInt);
                        toff2_tr = ton2_tr+toff1_tr(itVr)-ton1_tr(itVr);
                        idx = Target_tr(itVr).time{iInt}>ton2_tr;
                        iChR2_2p = iChR2fun_2p(p,X_tr(itVr,:),Target_tr(itVr).time{iInt},ton1_tr(itVR),toff1_tr(itVR),ton2_tr,toff2_tr);
                        iChR2_2p = iChR2_2p(idx);
                        iChR2_peakval = iChR2_2p(max(abs(iChR2_2p))==abs(iChR2_2p));

                        iChR2_measured = Target_tr(itVR).I{iInt};
                        iChR2_measured = iChR2_measured(idx);
                        iChR2_m_peakval = iChR2_measured(max(abs(iChR2_measured))==abs(iChR2_measured));

                        if iChR2_m_peakval~=0 && iChR2_peakval~=0
                            diffsq_mat = diffsq_mat+sum((iChR2_2p/iChR2_peakval-iChR2_measured/iChR2_m_peakval).^2);
                            nsamples = nsamples+length(iChR2_2p);
                            iChR2_in = iChR2_in + sum((iChR2_measured/iChR2_m_peakval).^2);
                        end
                    end
                end
            end
        case 'normalized+weighted'
            diffsq_mat = 0;
            nsamples = 0;
            iChR2_in = 0;
            for io=1:length(tV)

                iChR2_sp = iChR2fun_sp(p,X(io,:),Target_sp(io).time,ton1_sp(io),toff1_sp(io));
                iChR2_peakval = iChR2_sp(max(abs(iChR2_sp))==abs(iChR2_sp));

                iChR2_measured = Target_sp(io).I;
                iChR2_m_peakval = Target_sp(io).Ipeak;

                w = options.w(io,:);
                if iChR2_m_peakval~=0 && iChR2_peakval~=0
                    diffsq_mat = diffsq_mat+sum((w.*(iChR2_sp/iChR2_peakval-iChR2_measured/iChR2_m_peakval)).^2);
                    nsamples = nsamples+length(iChR2_sp);
                    iChR2_in = iChR2_in + sum((w.*iChR2_measured/iChR2_m_peakval).^2);
                end
            end
            if incl_tr_flag
                for itVr=1:length(X_tr)
                    for iInt=1:length(intervals)
                        ton2_tr = toff1_tr(itVR)+intervals(iInt);
                        toff2_tr = ton2_tr+toff1_tr(itVr)-ton1_tr(itVr);
                        idx = Target_tr(itVr).time{iInt}>ton2_tr;
                        iChR2_2p = iChR2fun_2p(p,X_tr(itVr,:),Target_tr(itVr).time{iInt},ton1_tr(itVR),toff1_tr(itVR),ton2_tr,toff2_tr);
                        iChR2_2p = iChR2_2p(idx);
                        iChR2_peakval = iChR2_2p(max(abs(iChR2_2p))==abs(iChR2_2p));

                        w_tr = options.w_tr(itVr,:);

                        iChR2_measured = Target_tr(itVR).I{iInt};
                        iChR2_measured = iChR2_measured(idx);
                        iChR2_m_peakval = iChR2_measured(max(abs(iChR2_measured))==abs(iChR2_measured));

                        if iChR2_m_peakval~=0 && iChR2_peakval~=0
                            diffsq_mat = diffsq_mat+sum((w_tr.*(iChR2_2p/iChR2_peakval-iChR2_measured/iChR2_m_peakval)).^2);
                            nsamples = nsamples+length(iChR2_2p);
                            iChR2_in = iChR2_in + sum((w_tr.*iChR2_measured/iChR2_m_peakval).^2);
                        end

                    end
                end
            end
        case 'zscore'
            diffsq_mat = 0;
            nsamples = 0;
            iChR2_in = 0;
            for io=1:length(tV)

                iChR2_sp = iChR2fun_sp(p,X(io,:),Target_sp(io).time,ton1_sp(io),toff1_sp(io));
                iChR2_measured = Target_sp(io).I;
                I_SD = Target_sp(io).I_SD;
                idx = I_SD~=0;

                iChR2_sp = iChR2_sp(idx);
                iChR2_measured = iChR2_measured(idx);
                I_SD = I_SD(idx);

                diffsq_mat = diffsq_mat+sum(((iChR2_sp-iChR2_measured)./I_SD).^2);
                nsamples = nsamples+length(iChR2_sp);
                iChR2_in = iChR2_in + sum((iChR2_measured./I_SD).^2);

            end
            if incl_tr_flag
                for itVr=1:length(X_tr)
                    for iInt=1:length(intervals)
                        ton2_tr = toff1_tr(itVR)+intervals(iInt);
                        toff2_tr = ton2_tr+toff1_tr(itVr)-ton1_tr(itVr);

                        I_SD = Target_tr(itVR).I_SD{iInt};
                        idx = Target_tr(itVr).time{iInt}>ton2_tr & I_SD~=0;

                        iChR2_2p = iChR2fun_2p(p,X_tr(itVr,:),Target_tr(itVr).time{iInt},ton1_tr(itVR),toff1_tr(itVR),ton2_tr,toff2_tr);
                        iChR2_2p = iChR2_2p(idx);
                        iChR2_measured = Target_tr(itVR).I{iInt};
                        iChR2_measured = iChR2_measured(idx);
                        I_SD = I_SD(idx);

                        diffsq_mat = diffsq_mat+sum(((iChR2_2p-iChR2_measured)./I_SD).^2);
                        nsamples = nsamples+length(iChR2_2p);
                        iChR2_in = iChR2_in + sum((iChR2_measured./I_SD).^2);
                    end
                end
            end
        case 'weighted'
            diffsq_mat = 0;
            nsamples = 0;
            iChR2_in = 0;
            for io=1:length(tV)

                iChR2_sp = iChR2fun_sp(p,X(io,:),Target_sp(io).time,ton1_sp(io),toff1_sp(io));
                iChR2_measured = Target_sp(io).I;

                w = options.w(io,:);

                diffsq_mat = diffsq_mat+sum((w.*(iChR2_sp-iChR2_measured)).^2);
                nsamples = nsamples+length(iChR2_sp);
                iChR2_in = iChR2_in + sum((w.*iChR2_measured).^2);

            end
            if incl_tr_flag
                for itVr=1:length(X_tr)
                    for iInt=1:length(intervals)
                        ton2_tr = toff1_tr(itVR)+intervals(iInt);
                        toff2_tr = ton2_tr+toff1_tr(itVr)-ton1_tr(itVr);

                        iChR2_2p = iChR2fun_2p(p,X_tr(itVr,:),Target_tr(itVr).time{iInt},ton1_tr(itVR),toff1_tr(itVR),ton2_tr,toff2_tr);
                        iChR2_measured = Target_tr(itVR).I{iInt};

                        w_tr = options.w_tr(itVr,:);

                        diffsq_mat = diffsq_mat+sum((w_tr.*(iChR2_2p-iChR2_measured)).^2);
                        nsamples = nsamples+length(iChR2_2p);
                        iChR2_in = iChR2_in + sum(w_tr.*(iChR2_measured).^2);
                    end
                end
            end
        otherwise
            error('false input')
    end

    RMS = sqrt(diffsq_mat/nsamples);
    Q2 = RMS/sqrt(iChR2_in);
    switch lower(options.obj)
        case 'rms'
            OBJ = RMS;
        case 'q2'
            OBJ = Q2;
        otherwise
            error ('false obj')
    end
catch
    OBJ = 1000;
end
end