function Target=gentTarget_surrogateData(funhandle,Vm_val,Il_val,OSpstart,OSpdsp,intervals,tend_sp,tend_2p,input,plot_flag)
% with this funciton we generate surrogate data to fit a 22HH model to. The
% function from which the target is extracted is given in funhandle as a
% string, Vm_val is a array with all the potentials of interest. Il val is
% a array with al light intensities of interest. OSpstart is the start of
% pulse in s, OSpdsp duration is duration of single pulse (DC*prp)
% Intervals is time in between pulses, OSd is total stimulation time, tend
% is the simulation end, input contains extra info needed for the funciton
% handle

Vm = repmat(Vm_val',1,length(Il_val)); Vm = Vm(:);
Il = repmat(Il_val,length(Vm_val),1); Il = Il(:);

if length(Vm)~=length(Il)
    error('lengths do not match')
end
if plot_flag
    figure(1)
    set(gcf,'position',[ -1717,105,1547,833]);
end

fh = str2func(funhandle);
% single pulse first
for i=1:length(Vm)
    Stimoptnames={'OSpstart','OSpd','OSdc','OSprf','OSipa','Tsim','V0'};
    Stimoptvalues=[OSpstart, OSpdsp, 1, tend_sp, Il(i), tend_sp,Vm(i)];
    for isov=1:length(Stimoptvalues)
        input.Stimopt.(Stimoptnames{isov})=Stimoptvalues(isov);
    end
    Out = fh(input);
    Features = Extract_feat(Out.time,Out.iChR2,OSpdsp,OSpstart,20,0);
    Target.SingleP(i).Vm = Vm(i);
    Target.SingleP(i).Il = Il(i);
    Target.SingleP(i).OSpstart = OSpstart;
    Target.SingleP(i).OSpd = OSpdsp;
    Target.SingleP(i).Ipeak = Features.Ipeak;
    Target.SingleP(i).Iss = Features.Iss;
    Target.SingleP(i).Iratio = Features.Iratio;
    Target.SingleP(i).TauOn = Features.TauOn;
    Target.SingleP(i).TauOff = Features.TauOff;
    Target.SingleP(i).TauInact = Features.TauInact;
    Target.SingleP(i).I = Out.iChR2;
    Target.SingleP(i).time = Out.time;
    Target.SingleP(i).nsamples = 1;

    if plot_flag
        ip = mod(i-1,length(Vm_val))+1;
        figure(1)
        subplot(ceil(length(Vm_val)/2),2,ip)
        hold on
        plot(Out.time,Out.iChR2,'DisplayName',sprintf('Il = %5.2f W/m�',Il(i)))
        hold off
        legend('show')
        xlabel('time (s)')
        ylabel('current');
        title(sprintf('Vm = %5.2f mV',Vm(i)))
    end

end

% double pulse
for i=1:length(Vm)

    Peak1 = zeros(1,length(intervals));
    Peak2 = zeros(1,length(intervals));
    Ratio = zeros(1,length(intervals));
    tend = zeros(1,length(intervals));
    time = cell(1,length(intervals));
    iChR2 = cell(1,length(intervals));

    for iInt = 1:length(intervals)

        OSprp = intervals(iInt)+OSpdsp;
        OSdc = OSpdsp/OSprp;
        Stimoptnames={'OSpstart','OSpd','OSdc','OSprf','OSipa','Tsim','V0'};
        Stimoptvalues=[OSpstart, 2*OSprp, OSdc, 1/OSprp, Il(i), tend_2p(iInt),Vm(i)];
        for isov=1:length(Stimoptvalues)
            input.Stimopt.(Stimoptnames{isov})=Stimoptvalues(isov);
        end
        Out = fh(input);
        time{iInt} = Out.time;
        iChR2{iInt} = Out.iChR2;
        tend(iInt) = Out.time(end);

        loc1=Out.time>=OSpstart & Out.time<=OSpstart+OSpdsp+1/2*intervals(iInt);
        loc2=Out.time>=OSpstart+OSprp & Out.time<=OSpstart+OSprp+OSpdsp+1/2*intervals(iInt);
        Peak1(iInt)=max(abs(Out.iChR2(loc1)));
        %Peak2(i)=max(abs(Out(i).iChR2(loc2)));
        Peaks2=findpeaks(abs(Out.iChR2(loc2)));
        Peak2(iInt)=Peaks2(1);
        Ratio(iInt)=Peak2(iInt)/Peak1(iInt);
    end

    xx=[0:1e-4:max(intervals)];
    RATIO=spline(intervals,Ratio,xx);
    [~,idx]=min(abs(RATIO-(1-exp(-1))));
    if any(Ratio<=(1-exp(-1))) && any(Ratio>=(1-exp(-1)))
        TauRecov=xx(idx(1));
    else
        TauRecov=nan;
    end

    Target.TauRecov(i).Vm = Vm(i);
    Target.TauRecov(i).Il = Il(i);
    Target.TauRecov(i).Intervals = intervals;
    Target.TauRecov(i).Pstart = OSpstart;
    Target.TauRecov(i).PD = OSpdsp;
    Target.TauRecov(i).nrpulses = 2;
    Target.TauRecov(i).TauRecov = TauRecov;
    Target.TauRecov(i).time = time;
    Target.TauRecov(i).I = iChR2;
    Target.TauRecov(i).tend = tend;
    Target.TauRecov(i).nsamples = 1;

    if plot_flag

        ip = floor((i-1)/length(Il_val))+1;
        figure(2)
        subplot(ceil(length(Vm_val)/2),2,ip)
        hold on
        scatter(Vm(i),TauRecov)
        hold off
        xlabel('Vm (mV)')
        ylabel('TauRecov(s)');
        title(sprintf('Il = %5.2f W/m�',Il(i)))
    end
end
end