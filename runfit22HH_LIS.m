function runfit22HH_LIS(filename,varargin)
% by running this function extra metadata is added in comparison with
% fit22HH. (LIS = Load Input and Settings). Aside from the input filename,
% three extra boolean inputs can be added (not required). First flag for
% saving plots generated by tool. Second, start profiler: analyze time
% spent in certain functions (profiler does not work when standalone applicaiton is created)
%Third and final, send mail when finished.

if nargin>3
    if varargin{3}
         %emailInfo contains a email adress, password and server
         % info can also be manually put in lines 88-90
         % change ofcourse according to own folders
        emailInfo = load('emailInfo.mat');
    end
end

if nargin>2
    if varargin{2}
        disp('profiler on')
        try
            mpiprofile on
        catch
            try
                profile on
            catch
                warning('no profiler')
                varargin{2} = 0;
            end
        end
    end
end

prefix = filename;
if contains(prefix,'.mat')
    prefix(end-3:end) = [];
end

suffix = datestr(now,'yymmddHH');
if contains(prefix,'_spsoBC_')
    suffix_extra = 'psoBC';
elseif contains(prefix,'_sms_')
    suffix_extra = 'ms';
elseif contains(prefix,'_sga_')
    suffix_extra = 'ga';
else
    suffix_extra = '';
end

% create folder where results are saved + add text file of progress (diary)
mkdir(sprintf('results%s%s',suffix,suffix_extra));
diary(sprintf('./results%s%s/run_%s_%s.txt',suffix,suffix_extra,prefix,suffix))

% parallel computing: run computation on all available cores according to
% 'local' profile
try
    c = parcluster('local');
    pool = parpool(c.NumWorkers);
end

try
Input = load(fullfile('./Inputs',filename));
catch
    Input = load(filename);
end

out = fit22HH(Input.Input,Input.settings);
save(sprintf('./results%s%s/run_%s_%s.mat',suffix,suffix_extra,prefix,suffix),'out','Input');
diary off

if nargin>1
    if varargin{1}
        disp('saving plots')
        %save plots
        number = get(gcf,'number');
        for i = 1:number
            figure(i);
            figname  = sprintf('./results%s%s/fig_%s_%s_%i.fig',suffix,suffix_extra,prefix,suffix,i);
            savefig(gcf,figname);
        end
    end
end
if nargin>2
    if varargin{2}
        disp('saving profiler')
        profile off
        profsave(profile('info'),sprintf('./results%s%s/profile_%s_%s',suffix,suffix_extra,prefix,suffix));
    end
end
if nargin>3
    if varargin{3}
        mail = emailInfo.mail;    % Replace with your email address
        password = emailInfo.password;       % Replace with your email password
        server = emailInfo.server;     % Replace with your SMTP server

        % Apply prefs and props
        props = java.lang.System.getProperties;
        props.setProperty('mail.smtp.port', '587');
        props.setProperty('mail.smtp.auth','true');
        props.setProperty('mail.smtp.starttls.enable','true');
        setpref('Internet','E_mail', mail);
        setpref('Internet','SMTP_Server', server);
        setpref('Internet','SMTP_Username', mail);
        setpref('Internet','SMTP_Password', password);
        % Send the email
        sendmail(mail,sprintf('finished_run_%s_%s',prefix,suffix),...
            sprintf('finished_run_%s_%s',prefix,suffix));
    end
end
pause(60)
exit
end