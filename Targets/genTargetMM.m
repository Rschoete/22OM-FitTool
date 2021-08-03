% generate target structure
% Target stucture should consist out ouf two fields: SingleP and TauRecov
% Idealy SingleP contains
%   Vm (mV)             Il(W/m�)
%   OSpstart (s)        OSpd(s)
%   Ipeak               Ipeak_SD
%   Iss                 Iss_SD
%   Iratio              Iratio_SD
%   TauOn(s)            TauOn_SD
%   TauOff(s)           TauOff_SD
%   TauInact(s)         TauInact_SD
%   I                   I_SD
%   time (s)            nsamples
% TauRecov contains:
%   Vm(mV)              Il(W/m�)
%   Intervals(s)(this is time inbetween pulses)
%   PD (s) pulse duration of single pulse (time on)
%   Pstart (s) onset of pulses
%   nrpulses
%   TauRecov(s)         TauRecov_SD
%   time (cell containing all time traces)
%   I cell containing all current traces
%   I_SD
%   tend                nsamples
clear all
close all
plot_flag = 1;
addpath(genpath('../../../Data_MM'))
addpath('../Functions')
ft = 1e-3; % input data in ms
%% IV files
filenames_IV = {'18628000','18628005','18628009'};%,'18628013'};
    [data, stimvals] = importfile_IV([filenames_IV{1},'.csv']);
    data_varNames = data.Properties.VariableNames;
    X0 = [0.001,0.001,0.001];
    LB = [0,0,0];
    UB = [0.05,0.1,0.5];
    nrruns = 100;
for ivar = 2:length(data_varNames)
    for iIV = 1:length(filenames_IV)
        [data, stimvals] = importfile_IV([filenames_IV{iIV},'.csv']);
        data_varNames = data.Properties.VariableNames;
        pause(0.1)
        OSPstart = 0.74376;
        OSPD = 0.5;
        figure
        plot(ft*data.(data_varNames{1}),data.(data_varNames{ivar}))
        Feature = Extract_feat(ft*data.(data_varNames{1}),data.(data_varNames{ivar}),OSPD,OSPstart,nrruns,plot_flag,...
            'X0',X0,'LB',LB,'UB',UB);
        collect.Imean(iIV,1) = Feature.Imean; collect.Iss(iIV,1) = Feature.Iss;
        collect.Ipeak(iIV,1) = Feature.Ipeak; collect.ttp(iIV,1) = Feature.ttp;
        collect.Iratio(iIV,1) = Feature.Iratio; collect.TauOn(iIV,1) = Feature.TauOn;
        collect.TauInact(iIV,1) = Feature.TauInact; collect.TauOff(iIV,1) = Feature.TauOff;
        collect.I(iIV,:) = data.(data_varNames{ivar});

    end
    if stimvals.V(ivar-1) == -60
        collect_60 = collect;
        place = ivar-1;
    end
    Target.SingleP(ivar-1).Vm = stimvals.V(ivar-1); Target.SingleP(ivar-1).OSpstart = OSPstart;
    Target.SingleP(ivar-1).Il = stimvals.Il(ivar-1);Target.SingleP(ivar-1).OSpd = OSPD;
    Target.SingleP(ivar-1).Ipeak = mean(collect.Ipeak); Target.SingleP(ivar-1).Iss = mean(collect.Iss);
    Target.SingleP(ivar-1).Ipeak_SD = std(collect.Ipeak); Target.SingleP(ivar-1).Iss_SD = std(collect.Iss);
    Target.SingleP(ivar-1).Imean = mean(collect.Imean); Target.SingleP(ivar-1).Iratio = mean(collect.Iratio);
    Target.SingleP(ivar-1).Imean_SD = std(collect.Imean); Target.SingleP(ivar-1).Iratio_SD = std(collect.Iratio);
    Target.SingleP(ivar-1).TauOn = mean(collect.TauOn); Target.SingleP(ivar-1).TauOff = mean(collect.TauOff);
    Target.SingleP(ivar-1).TauOn_SD = std(collect.TauOn); Target.SingleP(ivar-1).TauOff_SD = std(collect.TauOff);
    Target.SingleP(ivar-1).TauInact = nanmean(collect.TauInact); Target.SingleP(ivar-1).ttp = mean(collect.ttp);
    Target.SingleP(ivar-1).TauInact_SD = nanstd(collect.TauInact); Target.SingleP(ivar-1).ttp_SD = std(collect.ttp);
    Target.SingleP(ivar-1).I = mean(collect.I,1); Target.SingleP(ivar-1).I_SD = std(collect.I,[],1);
    Target.SingleP(ivar-1).time = ft*data.(data_varNames{1})'; Target.SingleP(ivar-1).nsamples = length(filenames_IV);
    collect = struct();
    pause(0.01)
end
NIV = length(Target.SingleP);
%% LightTitration
filenames_Il = {'18524054', '18524058'};%, '18606021', '18606025'};
Fs = 1e4;
Niir = 30;
Fst = 200;
plotnr = 8;
ycol = [];
[data, stimvals] = importfile_Il([filenames_Il{1},'.csv']);
data_varNames = data.Properties.VariableNames;
for ivar = 2:length(data_varNames)
    for iIl = 1:length(filenames_Il)
        [data, stimvals] = importfile_Il([filenames_Il{iIl},'.csv']);
        data_varNames = data.Properties.VariableNames;
        pause(0.1)
        OSPD = 2;
        OSPstart = 0.893;
        figure
        plot(ft*data.(data_varNames{1}),data.(data_varNames{ivar}))
        Feature = Extract_feat(ft*data.(data_varNames{1}),data.(data_varNames{ivar}),OSPD,OSPstart,nrruns,plot_flag,...
            'X0',X0,'LB',LB,'UB',UB);
        collect.Imean(iIl,1) = Feature.Imean; collect.Iss(iIl,1) = Feature.Iss;
        collect.Ipeak(iIl,1) = Feature.Ipeak; collect.ttp(iIl,1) = Feature.ttp;
        collect.Iratio(iIl,1) = Feature.Iratio; collect.TauOn(iIl,1) = Feature.TauOn;
        collect.TauInact(iIl,1) = Feature.TauInact; collect.TauOff(iIl,1) = Feature.TauOff;
        collect.I(iIl,:) = data.(data_varNames{ivar});
    end
    if abs(stimvals.Il(ivar-1)-3.734436e3)<1e-6 && 0
        fnms = fieldnames(collect);
        fnms(end) = [];
        cell1 = struct2cell(collect); cell1(end) = [];
        cell2 = struct2cell(collect_60); cell2(end) = [];
        collect_all = cell2struct(cellfun(@vertcat,cell1,cell2,'uni',0),fnms,1);
        collect_intm = collect;
        collect = collect_all; collect.I = collect_intm.I;
        idx_pos = NIV+ivar-1;
    else
        idx_pos = NIV+ivar-1;
    end
    Target.SingleP(idx_pos).Vm = stimvals.V(ivar-1); Target.SingleP(idx_pos).OSpstart = OSPstart;
    Target.SingleP(idx_pos).Il = stimvals.Il(ivar-1);Target.SingleP(idx_pos).OSpd = OSPD;
    Target.SingleP(idx_pos).Ipeak = mean(collect.Ipeak); Target.SingleP(idx_pos).Iss = mean(collect.Iss);
    Target.SingleP(idx_pos).Ipeak_SD = std(collect.Ipeak); Target.SingleP(idx_pos).Iss_SD = std(collect.Iss);
    Target.SingleP(idx_pos).Imean = mean(collect.Imean); Target.SingleP(idx_pos).Iratio = mean(collect.Iratio);
    Target.SingleP(idx_pos).Imean_SD = std(collect.Imean); Target.SingleP(idx_pos).Iratio_SD = std(collect.Iratio);
    Target.SingleP(idx_pos).TauOn = mean(collect.TauOn); Target.SingleP(idx_pos).TauOff = mean(collect.TauOff);
    Target.SingleP(idx_pos).TauOn_SD = std(collect.TauOn); Target.SingleP(idx_pos).TauOff_SD = std(collect.TauOff);
    Target.SingleP(idx_pos).TauInact = mean(collect.TauInact); Target.SingleP(idx_pos).ttp = mean(collect.ttp);
    Target.SingleP(idx_pos).TauInact_SD = std(collect.TauInact); Target.SingleP(idx_pos).ttp_SD = std(collect.ttp);
    Target.SingleP(idx_pos).I = mean(collect.I,1); Target.SingleP(idx_pos).I_SD = std(collect.I,[],1);
    Target.SingleP(idx_pos).time = ft*data.(data_varNames{1})'; Target.SingleP(idx_pos).nsamples = length(collect.Imean);
    collect = struct();
    pause(0.01)
end


%% Recov
filenames_recov = {'18620032', '18801001', '18808002', '18808014'};
OSPstart = 1.0905;
OSPD = 1;
[data, stimvals] = importfile_recov([filenames_recov{1},'.csv']);
data_varNames = data.Properties.VariableNames;
for irecov = 1:length(filenames_recov)
    [data, stimvals] = importfile_recov([filenames_recov{irecov},'.csv']);
    data_varNames = data.Properties.VariableNames;
    pause(0.01)
    for ivar = 2:length(data_varNames)
        figure
        plot(ft*data.(data_varNames{1}),data.(data_varNames{ivar}))
        loc1=ft*data.(data_varNames{1})>=OSPstart & ft*data.(data_varNames{1})<=OSPstart+OSPD+1/2*stimvals.dt(ivar-1);
        loc2=ft*data.(data_varNames{1})>=OSPstart+OSPD+stimvals.dt(ivar-1) & ...
            ft*data.(data_varNames{1})<=OSPstart+2*OSPD+3/2*stimvals.dt(ivar-1);
        Peak1(ivar-1)=max(abs(data.(data_varNames{ivar})(loc1)));
        %Peak2(i)=max(abs(Out(i).iChR2(loc2)));
        Peaks2=findpeaks(abs(data.(data_varNames{ivar})(loc2)));
        Peak2(ivar-1)=max(Peaks2);
        Ratio(ivar-1)=Peak2(ivar-1)/Peak1(ivar-1);
        pause(0.01)
        idx_toretain = ft*data.(data_varNames{1})'<=10+OSPstart+OSPD+stimvals.dt(ivar-1); %retain 10s after start second pulse
        collect.I(irecov).cell{ivar-1} = data.(data_varNames{ivar})(idx_toretain)';
        collect.time(irecov).cell{ivar-1} = ft*data.(data_varNames{1})(idx_toretain)';
        collect.tend(irecov).tend(ivar-1) = collect.time(irecov).cell{ivar-1}(end);
    end
    xx=[0:1e-4:max(stimvals.dt)];
    RATIO=spline(stimvals.dt(1:end),Ratio,xx);
    [~,idx]=min(abs(RATIO-(1-exp(-1))));
    if any(Ratio<=(1-exp(-1))) && any(Ratio>=(1-exp(-1)))
        TauRecov=xx(idx(1));
    else
        TauRecov=nan;
    end
    collect.TauRecov(irecov) = TauRecov;
    clearvars('data','Ratio','Peak1','Peak2')
end
Target.TauRecov.Vm = stimvals.V(1);
Target.TauRecov.Il = stimvals.Il(1);
Target.TauRecov.Intervals = stimvals.dt;
Target.TauRecov.PD = OSPD;
Target.TauRecov.Pstart = OSPstart;
Target.TauRecov.nrpulses = 2;
Target.TauRecov.TauRecov = mean(collect.TauRecov);
Target.TauRecov.TauRecov_SD = std(collect.TauRecov);
Target.TauRecov.time = collect.time(1).cell;
Target.TauRecov.tend = collect.tend(1).tend;
Target.TauRecov.nsamples = length(filenames_recov);
for ivar = 1:length(data_varNames)-1
    Iall = [];
    for irecov = 1:length(filenames_recov)
    Iall = vertcat(Iall,collect.I(irecov).cell{ivar});
    end
    Isave{ivar} = mean(Iall,1);
    Isave_SD{ivar} = std(Iall,[],1);
end
Target.TauRecov.I = Isave;
Target.TauRecov.I_SD = Isave_SD;

if input('save (1/0)')

    save(['TargetMM',datestr(now,'yymmdd'),'.mat'],'Target')
end