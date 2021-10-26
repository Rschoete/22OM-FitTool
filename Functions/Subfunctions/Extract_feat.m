function [Feature] = Extract_feat(time,I,OSpd,OSpstart,nrruns,PLOT,varargin)
% this function extracts features necessary for fitting opsins
%Ipeak,Iss,Iratio,Tauon,Tauoff,TauInact
method = 'both';
ssReached_flag = 1; %is steady-state reached at end of pulse? if False steady state value is extracted via curve fit
segTauInact_startpoint = 0.005; %shift wrt peak current indicating start of segment Tau inact
global fignr
if isempty(fignr)
    fignr = get(gcf,'number')+1;
    idx_sp = 0;
else
    idx_sp = 1;
end

if license('test','Parallel_Computing_Toolbox') %|| license('test','Distrib_computing_Toolbox')
    Parallel_flag = 'always';
else
    Parallel_flag = 0;
end

X0=0.01; LB=0; UB=1.5;
if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'X0'))
            X0 = varargin{find(strcmpi(varargin,'X0'))+1};
        end
        if any(strcmpi(varargin,'LB'))
            LB = varargin{find(strcmpi(varargin,'LB'))+1};
        end
        if any(strcmpi(varargin,'UB'))
            UB = varargin{find(strcmpi(varargin,'UB'))+1};
        end
        if any(strcmpi(varargin,'method'))
            method = varargin{find(strcmpi(varargin,'method'))+1};
        end
        if any(strcmpi(varargin,'ssReached_flag'))
            % is steady-state (ss) value reached at the end of pulse?
            ssReached_flag = varargin{find(strcmpi(varargin,'ssReached_flag'))+1};
        end
        if any(strcmpi(varargin,'segTauInact_startpoint'))
            % is steady-state (ss) value reached at the end of pulse?
            segTauInact_startpoint = varargin{find(strcmpi(varargin,'segTauInact_startpoint'))+1};
        end
    end
end
%initial parameters
EXP = @(b,t) b(1).*exp(-t./b(2))+b(3);
Logistic = @(b,t) b(1)./(1+exp(t./b(2)))+b(3);
OBJ = {EXP,EXP,EXP};
Names = {'Mono-exp','Mono-exp','Mono-exp'};
TauNames = {'\tau_{on}','\tau_{off}','\tau_{inact}'};
%h=waitbar(0,'feature extraction');
opt=optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'MaxFunctionEvaluations',600e3,'MaxIterations',600e3);
STOP=3;
%Feature.Imean = trapz(time(time>=OSpstart),I(time>=OSpstart))/(time(end)-OSpstart);
%Segmentation for curve fitting
[~,loc_peaks] = max(abs(I(time>=OSpstart & time<=OSpstart+OSpd))); %find peak (during stimulation)
loc_peaks = loc_peaks+sum(time<OSpstart); %adjust to right position
Feature.Ipeak = I(loc_peaks(1)); % place peak in structure current
Iss_pulseEnd = mean(I(time>=OSpstart+0.9*OSpd & time<=OSpstart+0.95*OSpd)); % get steady state current candidate
Feature.ttp = time(loc_peaks(1))-OSpstart;
if isnan(Iss_pulseEnd )
    Iss_pulseEnd  = I(time==max(time(time<=OSpstart+OSpd))); %if previous method couldnt define SS now take current closest but smaller than end of pulse
end
if abs(Feature.Ipeak)<=abs(Iss_pulseEnd) %Peak must always be larger than SS
    STOP = 2; % no inactivation curve
    Feature.Ipeak = Iss_pulseEnd;
    %loc_peaks=find(time==max(time(time<=OSpstart+OSpd)));
    loc_peaks = find(time>OSpstart & time<=(OSpstart+OSpd) & I==Iss_pulseEnd);
end
if ssReached_flag
    Feature.Iss = Iss_pulseEnd;
    Feature.Iratio=Feature.Iss/Feature.Ipeak;
end
%Tau on segment
%SegTauon=find(time>=OSpstart+0.2*Feature.ttp & time<=time(loc_peaks(1))); % identify indices for segment Tauon
SegTauonstart = find(abs(I)<=abs(Feature.Ipeak-0.75*(Feature.Ipeak-I(abs(time-OSpstart)==min(abs(time-OSpstart)))))); 
%this method eliminates inclusion of photon reaction time see paper
SegTauonstart = SegTauonstart(SegTauonstart<loc_peaks(1));
SegTauon=find(time>=time(SegTauonstart(end)) & time<=time(loc_peaks(1))); % identify indices for segment Tauon
correct_flag = 1; % we correct becasue photon reaction time is not included
if length(SegTauon)<3
% if insufficient timepoints increase Tauon segment
    SegTauon=find(time>=OSpstart & time<=time(loc_peaks(1)));
    correct_flag = 0;
end
if length(SegTauon)<3
    SegTauon=find(time>=0.1*OSpstart & time<=time(loc_peaks(1)));
    correct_flag = 0;
end
tTauon=(time(SegTauon)-time(SegTauon(1))); ITauon=I(SegTauon);% time axis for tau on



%Tau off segment
SegTauoff=find(time>=(OSpstart+OSpd) & time<=OSpstart+1.2*OSpd);
if length(SegTauoff)<3
    SegTauoff=find(time>=(OSpstart+OSpd) & abs(I)>=0);
end
tTauoff=(time(SegTauoff)-time(SegTauoff(1))); ITauoff=I(SegTauoff);

%Tau inact segment
if time(loc_peaks(1))<=OSpstart+0.90*OSpd && STOP==3
    try
        SegTauinact=find(time>=(time(loc_peaks(1))+segTauInact_startpoint) & time<=OSpstart+0.99*OSpd);
        tTauinact=(time(SegTauinact)-time(SegTauinact(1))); ITauinact=I(SegTauinact);
        STOP=3;
        if length(tTauinact)<3
            try
                SegTauinact=find(time>=(time(loc_peaks(1))+segTauInact_startpoint) & time<=(OSpstart.OSpd));
                tTauinact=(time(SegTauinact)-time(SegTauinact(1))); ITauinact=I(SegTauinact);
            catch
                tTauinact=[]; ITauinact=[];
                STOP=2;
                Feature.TauInact=NaN;
            end
        end
    catch
        tTauinact=[]; ITauinact=[];
        STOP=2;
        Feature.TauInact=NaN;
    end
else
    tTauinact=[]; ITauinact=[];
    STOP=2;
    Feature.TauInact=NaN;
end
%Gather in list
tTau={tTauon,tTauoff,tTauinact}; ITau={ITauon,ITauoff,ITauinact};


switch lower(method)
    case {'both','fit'}
        b1=[-(ITauon(end)-ITauon(1)),Iss_pulseEnd];
        b3=[ITauon(end),0];
        % add seperatly because Itauinact can be empty
        if STOP==3 && ssReached_flag
            b1(3)=-(Feature.Iss-ITauinact(1));
            b3(3)=Feature.Iss;
        end
        if length(X0)==1
            X0 = repmat(X0,1,STOP);
        end
        if length(LB)==1
            LB = repmat(LB,1,STOP);
        end
        if length(UB)==1
            UB = repmat(UB,1,STOP);
        end

        for i=1:STOP %loop over different tau (all same OBJ)
            %waitbar(i/STOP,h,TauNames{i});

            if i==3 && ~ssReached_flag
                LB_intm = [0,LB(i),0];
                UB_intm = [-10*(Iss_pulseEnd-ITauinact(1)),UB(i),Iss_pulseEnd*10];
                UB_intm2 = UB_intm;
                UB_intm(UB_intm2<LB_intm)=LB_intm(UB_intm2<LB_intm);
                LB_intm(UB_intm2<LB_intm)=UB_intm2(UB_intm2<LB_intm);
                objfun=@(p,t) OBJ{i}([p(1),p(2),p(3)],t);
                problem=createOptimProblem('lsqcurvefit',...
                    'objective',objfun,...
                    'xdata',tTau{i},'ydata',ITau{i},...
                    'x0',[-(Iss_pulseEnd/2-ITauinact(1)),X0(i),Iss_pulseEnd/2],...
                    'lb',LB_intm,...
                    'ub',UB_intm, 'options',opt);
                ms = MultiStart('Display','off','UseParallel',Parallel_flag,'FunctionTolerance',1e-12,'XTolerance',1e-12);
                [b,fval,~,~,solutions]= run(ms,problem,nrruns);
                Feature.Iss = b(3);
                Feature.Iratio=Feature.Iss/Feature.Ipeak;
                b = b(2);

            else
                objfun=@(p,t) OBJ{i}([b1(i),p,b3(i)],t);
                    problem=createOptimProblem('lsqcurvefit',...
                    'objective',objfun,...
                    'xdata',tTau{i},'ydata',ITau{i},...
                    'x0',X0(i),...
                    'lb',LB(i),...
                    'ub',UB(i), 'options',opt);
                ms = MultiStart('Display','off','UseParallel',Parallel_flag,'FunctionTolerance',1e-12,'XTolerance',1e-12);
                [b,fval,~,~,solutions]= run(ms,problem,nrruns);
            end
            %pause(0.01)
            if i==1
                if correct_flag
                    b = b*(1-log((I(SegTauon(end))-I(SegTauon(1)))/I(SegTauon(end))));
                end
                Feature.TauOn=b;

            elseif i==3
                Feature.TauInact=b;
            else
                Feature.TauOff=b;
            end
            Feature.TauSol{i}=b;
            Feature.TauSSres(i)=fval;
            Feature.TauR2(i)=1-fval/sum((ITau{i}-mean(ITau{i})).^2);
            if PLOT
                Colors=[0, 0.4470, 0.7410;
                    0.8500, 0.3250, 0.0980;
                    0.9290, 0.6940, 0.1250;
                    0.4940, 0.1840, 0.5560;
                    0.4660, 0.6740, 0.1880;
                    0.3010, 0.7450, 0.9330;
                    0.6350, 0.0780, 0.1840
                    0,0,1;
                    1,0,0;
                    0,1,0;
                    1,1,0;
                    1,0,1];
                fig(fignr)
                subplot(2,2,i)
                plot(tTau{i},ITau{i},'.','color',Colors(2,:),'Displayname','Data')
                hold on
                if ~(i==3 && ~ssReached_flag)
                    % wont work in case ssReached_flag false because objfun is different then!
                    plot(tTau{i},objfun(solutions(1,1).X,tTau{i}),'color',Colors(1,:),'DisplayName',[Names{i},', \tau = ',num2str(solutions(1,1).X(1))])
                end
                title(TauNames{i});

                drawnow
            end
        end
        if strcmpi(method,'both')
            IOI = (1-exp(-1))*Feature.Ipeak;
            Feature.TauOn2 = time(SegTauon(min(abs(I(SegTauon)-IOI))==abs(I(SegTauon)-IOI)))-OSpstart;
            IOI = exp(-1)*Iss_pulseEnd;
            Feature.TauOff2 = time(SegTauoff(min(abs(I(SegTauoff)-IOI))==abs(I(SegTauoff)-IOI)))-time(SegTauoff(1));
            if STOP>2
                IOI = exp(-1)*(Feature.Ipeak-Feature.Iss)+Feature.Iss;
                Feature.TauInact2 = time(SegTauinact(min(abs(I(SegTauinact)-IOI))==abs(I(SegTauinact)-IOI)))-time(loc_peaks(1));
            end
        end
    case 'percentage'
        IOI = (1-exp(-1))*Feature.Ipeak;
        Feature.TauOn = time(SegTauon(min(abs(I(SegTauon)-IOI))==abs(I(SegTauon)-IOI)))-OSpstart;
        IOI = exp(-1)*Feature.Iss;
        Feature.TauOff = time(SegTauoff(min(abs(I(SegTauoff)-IOI))==abs(I(SegTauoff)-IOI)))-time(SegTauoff(1));
        if STOP>2
            if ssReached_flag
                IOI = exp(-1)*(Feature.Ipeak-Feature.Iss)+Feature.Iss;
                Feature.TauInact = time(SegTauinact(min(abs(I(SegTauinact)-IOI))==abs(I(SegTauinact)-IOI)))-time(loc_peaks(1));
            else
                Feature.Iss = Iss_pulseEnd;
                Feature.Iratio = Feature.Iss/Feature.Ipeak;
                IOI = exp(-1)*(Feature.Ipeak-Feature.Iss)+Feature.Iss;
                Feature.TauInact = time(SegTauinact(min(abs(I(SegTauinact)-IOI))==abs(I(SegTauinact)-IOI)))-time(loc_peaks(1));
                warning('steady state not reached => percentage might not give correct results')
            end
        end

    otherwise
        error ('wrong method')
end
if STOP==2 && ~ssReached_flag
    Feature.TauInact=NaN;
    Feature.Iss = NaN;
    Feature.Iratio = NaN;
end
%close(h)
end

