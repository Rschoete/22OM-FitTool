function [fignr1,fignr2] = plotFeatures(Target,Il_val,Vm_val,fTnames,fInames,Colors,errorbar_flag,...
    Target_tr,Il_valtr,Vm_valtr,plotTauRecov_flag,varargin)
%plot features versus potetnial V and light intensity Il

figpos_flag = 0;
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

for ifIn=1:length(fInames)
    % plot fInames wrt V
    subplot(length(fInames),2,(ifIn-1)*2+1);
    for iIl=1:length(Il_val)
        yvals = [Target([Target.Il]==Il_val(iIl)).(fInames{ifIn})];
        xvals = [Target([Target.Il]==Il_val(iIl)).Vm];
        idx = ~isnan(yvals);
        xvals = xvals(idx);
        yvals = yvals(idx);

        if length(yvals)>1
            if errorbar_flag
                yvals_SD = [Target([Target.Il]==Il_val(iIl)).([fInames{ifIn},'_SD'])];
                yvals_SD = yvals_SD(idx);
                errorbar(xvals,yvals,...
                    yvals_SD,'*','color',Colors(iIl),'DisplayName',sprintf('%5.2e',Il_val(iIl)));
            else
                plot(xvals,yvals,'*','color',Colors(iIl),'DisplayName',sprintf('%5.2e',Il_val(iIl)))
            end
        end
        hold on
    end
    hold off
    xlabel('Vm (mV)');
    title([fInames{ifIn}(1),'_{',fInames{ifIn}(2:end),'}']);
    legend('show')

    % plot Ifnmaes wrt Il
    subplot(length(fInames),2,ifIn*2);
    for iVm=1:length(Vm_val)
        yvals = [Target([Target.Vm]==Vm_val(iVm)).(fInames{ifIn})];
        xvals = [Target([Target.Vm]==Vm_val(iVm)).Il];
        idx = ~isnan(yvals);
        yvals = yvals(idx);
        xvals = xvals(idx);

        if length(yvals)>1
            if errorbar_flag
                yvals_SD = [Target([Target.Vm]==Vm_val(iVm)).([fInames{ifIn},'_SD'])];
                yvals_SD = yvals_SD(idx);
                errorbar(xvals,yvals,...
                    yvals_SD,'*',...
                    'color',Colors(iVm),'DisplayName',sprintf('%5.2e',Vm_val(iVm)));
            else
                plot(xvals,yvals,'*',...
                    'color',Colors(iVm),'DisplayName',sprintf('%5.2e',Vm_val(iVm)));
            end
        end
        hold on
    end
    set(gca,'Xscale','log')
    title([fInames{ifIn}(1),'_{',fInames{ifIn}(2:end),'}']);
    xlabel('Il (W/m�)');
    legend('show')
end

fignr2 = fignr1+1;
figure(fignr2);
if figpos_flag
    set(gcf,'position',figposition)
end

for ifTn=1:length(fTnames)
    % plot fTnames wrt V
    subplot(length(fTnames)+plotTauRecov_flag,2,(ifTn-1)*2+1);
    for iIl=1:length(Il_val)
        yvals = [Target([Target.Il]==Il_val(iIl)).(fTnames{ifTn})];
        xvals = [Target([Target.Il]==Il_val(iIl)).Vm];
        idx = ~isnan(yvals);
        yvals = yvals(idx);
        xvals = xvals(idx);

        if length(yvals)>1
            if errorbar_flag
                yvals_SD = [Target([Target.Il]==Il_val(iIl)).([fTnames{ifTn},'_SD'])];
                yvals_SD = yvals_SD(idx);
                errorbar(xvals,yvals,...
                    yvals_SD,'*','color',Colors(iIl),'DisplayName',sprintf('%5.2e',Il_val(iIl)));
            else
                plot(xvals,yvals,'*','color',Colors(iIl),'DisplayName',sprintf('%5.2e',Il_val(iIl)))
            end
        end
        hold on
    end
    hold off
    xlabel('Vm (mV)');
    title([fTnames{ifTn}(1:3),'_{',fTnames{ifTn}(4:end),'}']);
    legend('show')

    % plot Tfnmaes wrt Il
     subplot(length(fTnames)+plotTauRecov_flag,2,ifTn*2);
    for iVm=1:length(Vm_val)
        yvals = [Target([Target.Vm]==Vm_val(iVm)).(fTnames{ifTn})];
        xvals = [Target([Target.Vm]==Vm_val(iVm)).Il];
        idx = ~isnan(yvals);
        yvals = yvals(idx);
        xvals = xvals(idx);

        if length(yvals)>1
            if errorbar_flag
                yvals_SD = [Target([Target.Vm]==Vm_val(iVm)).([fTnames{ifTn},'_SD'])];
                yvals_SD = yvals_SD(idx);
                errorbar(xvals,yvals,...
                    yvals_SD,'*',...
                    'color',Colors(iVm),'DisplayName',sprintf('%5.2e',Vm_val(iVm)));
            else
                plot(xvals,yvals,'*',...
                    'color',Colors(iVm),'DisplayName',sprintf('%5.2e',Vm_val(iVm)));
            end
        end
        hold on
    end
    set(gca,'Xscale','log')
    title([fTnames{ifTn}(1:3),'_{',fTnames{ifTn}(4:end),'}']);
    xlabel('Il (W/m�)');
    legend('show')

end

if plotTauRecov_flag
   %plot TauRecov wrt V
    subplot(length(fTnames)+plotTauRecov_flag,2,(ifTn)*2+1);
    for iIl=1:length(Il_valtr)
        yvals = [Target_tr([Target_tr.Il]==Il_valtr(iIl)).('TauRecov')];
        xvals = [Target_tr([Target_tr.Il]==Il_valtr(iIl)).Vm];
        idx = ~isnan(yvals);
        yvals = yvals(idx);
        xvals = xvals(idx);
        if length(yvals)>=1
            if errorbar_flag && isfield(Target_tr,'TauRecov_SD')
                yvals_SD = [Target_tr([Target_tr.Il]==Il_valtr(iIl)).('TauRecov_SD')];
                yvals_SD = yvals_SD(idx);
                errorbar(xvals,yvals,...
                    yvals_SD,'*','color',Colors(iIl),'DisplayName',sprintf('%5.2e',Il_valtr(iIl)));
            else
                plot(xvals,yvals,'*','color',Colors(iIl),'DisplayName',sprintf('%5.2e',Il_valtr(iIl)))
            end
        end
        hold on
    end
    hold off
    xlabel('Vm (mV)');
    title('Tau_{Recov}');
    legend('show')

    % plot TauRecov wrt Il
    subplot(length(fTnames)+plotTauRecov_flag,2,ifTn*2+2);
    for iVm=1:length(Vm_valtr)
        yvals = [Target_tr([Target_tr.Vm]==Vm_valtr(iVm)).('TauRecov')];
        xvals = [Target_tr([Target_tr.Vm]==Vm_valtr(iVm)).Il];
        idx = ~isnan(yvals);
        yvals = yvals(idx);
        xvals = xvals(idx);

        if length(yvals)>=1
            if errorbar_flag && isfield(Target_tr,'TauRecov_SD')
                yvals_SD = [Target_tr([Target_tr.Vm]==Vm_valtr(iVm)).('TauRecov_SD')];
                yvals_SD = yvals_SD(idx);
                errorbar(xvals,yvals,...
                    yvals_SD,'*',...
                    'color',Colors(iVm),'DisplayName',sprintf('%5.2e',Vm_valtr(iVm)));
            else

                plot(xvals,yvals,'*',...
                    'color',Colors(iVm),'DisplayName',sprintf('%5.2e',Vm_valtr(iVm)));
            end
        end
        hold on
    end
    set(gca,'Xscale','log')
    title('Tau_{Reocv}');
    xlabel('Il (W/m�)');
    legend('show')
end
end