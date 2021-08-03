function plotTau(xdata,ydata,Vm_val,IlTau_val,errorbar_flag,varargin)
% plot training data of TauO and TauDA

figpos_flag = 0;
fns = {'TauO','TauDA'};

fignr1 = get(gcf,'number')+1;
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

figure(fignr1);
if figpos_flag
    set(gcf,'position',figposition)
end

for ifns =1:length(fns)
    subplot(length(fns),2,(ifns-1)*2+1)

    for iIl = 1:length(IlTau_val)
        idx = xdata.(fns{ifns})(:,2)==IlTau_val(iIl);
        yvals = ydata.(fns{ifns})(idx);
        xvals = xdata.(fns{ifns})(idx,1);
        if length(yvals)>1
            if errorbar_flag
                errorbar(xvals,yvals,ydata.([fns{ifns},'_SD'])(idx),'*','DisplayName',sprintf('%5.2e',IlTau_val(iIl)));
            else
                plot(xvals,yvals,'*','DisplayName',sprintf('%5.2e',IlTau_val(iIl)));
            end
            hold on
        end
    end

    hold off
    title([fns{ifns}(1:3),'_{',fns{ifns}(4:end),'}'])
    xlabel('Vm [mV]')
    legend show

    subplot(length(fns),2,ifns*2)

    for iVm = 1:length(Vm_val)
        idx = xdata.(fns{ifns})(:,1)==Vm_val(iVm);
        yvals = ydata.(fns{ifns})(idx);
        xvals = xdata.(fns{ifns})(idx,2);
        if length(yvals)>1
            if errorbar_flag
                errorbar(xvals,yvals,ydata.([fns{ifns},'_SD'])(idx),'*','DisplayName',sprintf('%5.2e',Vm_val(iVm)));
            else
                plot(xvals,yvals,'*','DisplayName',sprintf('%5.2e',Vm_val(iVm)));
            end
            hold on
        end
    end

    hold off
    title([fns{ifns}(1:3),'_{',fns{ifns}(4:end),'}'])
    xlabel('Il [mW/mm�]')
    legend show
end
end