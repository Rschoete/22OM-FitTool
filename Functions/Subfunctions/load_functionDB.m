function [fun,OinfIfun_names,TauOIfun_names,DAinfIfun_names,...
    TauDAIfun_names,TauOVfun_names,TauDAVfun_names,Grectfun_names, combTauO,combTauDA] = ...
    load_functionDB(varargin)
% create fitDatabase

set = 'standard';
newfunctions_flag = 0;
newnames_flag = 0;

% Change standard input parameters
if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'set'))
            set = varargin{find(strcmpi(varargin,'set'))+1};
        end
        if any(strcmpi(varargin,'newfun'))
            newfun = varargin{find(strcmpi(varargin,'newfun'))+1};
            newfunctions_flag = 1;
        end
        if ~isempty(intersect({varargin{1:2:end}},{'OinfI','TauOI','DAinfI','TauDAI','TauOV','TauDAV','Grect','combTauO','combTauDa'}))
            newnames_flag = 1;
        end
    end
end

% Il=0 => 0 Il=inf  => 1  bend at 10^p(1) slope at bend =  1/(4*p(2))/dec
fun.symLogisticslog.fun=@(p,I) (1./(1+exp(p(1)./p(2)).*I.^(-1./(p(2)*log(10)))));
fun.symLogisticslog.LB=[-10 0];
fun.symLogisticslog.X0=[1,1];
fun.symLogisticslog.UB=[10,20];
fun.symLogisticslog.nrpara=2;
fun.symLogisticslog.var=2;

% Il=0 => 0 Il=inf  => p(3) bend at 10^p(1) slope at bend =  p(3)/(4*p(2))/dec
fun.symLogisticsloga.fun=@(p,I) (p(3)./(1+exp(p(1)./p(2)).*I.^(-1./(p(2)*log(10)))));
fun.symLogisticsloga.LB=[-10,0,0];
fun.symLogisticsloga.X0=[1,1,1];
fun.symLogisticsloga.UB=[10,20,100];
fun.symLogisticsloga.nrpara=3;
fun.symLogisticsloga.var=2;

% Il=0 => 1 Il=inf  =>1-p(3) bend at 10^p(1) slope at bend =  -p(3)/(4*p(2))/dec
fun.symLogisticslogmina.fun=@(p,I) 1-(p(3)./(1+exp(p(1)/p(2)).*I.^(-1./(p(2)*log(10)))));
fun.symLogisticslogmina.LB=[-10,0,0];
fun.symLogisticslogmina.X0=[1,1,0.1];
fun.symLogisticslogmina.UB=[10,20,1];
fun.symLogisticslogmina.nrpara=3;
fun.symLogisticslogmina.var=2;

% Il=0 => 0 Il=inf  => inf init slope p2/p1
fun.logarithm.fun=@(p,I) p(2).*log(1+I./10^(p(1)));
fun.logarithm.nrpara=2;
fun.logarithm.LB=[-10 0];
fun.logarithm.X0=[1 1];
fun.logarithm.UB=[10,1e2];
fun.logarithm.var=2;

% Il=0 => p(3) Il=inf  => inf init slope p2/p1
fun.logarithm3.fun=@(p,I) p(2).*log(1+I./10^(p(1)))+p(3);
fun.logarithm3.nrpara=3;
fun.logarithm3.LB=[-10 0 0];
fun.logarithm3.X0=[1 1 0];
fun.logarithm3.UB=[10 100 10];
fun.logarithm3.var=2;

% Il = 0 ==> 0 Il=inf ==> inf
fun.polynomial.fun=@(p,I) p(2).*I.^p(1);
fun.polynomial.nrpara=2;
fun.polynomial.LB=[0 0];
fun.polynomial.X0=[1 1];
fun.polynomial.UB=[5 10];
fun.polynomial.var=2;

% Il=0 => 0 Il==inf ==> inf init slope = 0
fun.logpolynomial.fun=@(p,I) p(2).*log(1+I).^p(1);
fun.logpolynomial.nrpara=2;
fun.logpolynomial.LB=[0 0];
fun.logpolynomial.X0=[1 1];
fun.logpolynomial.UB=[10 10];
fun.logpolynomial.var=2;

% Il=0 => p(3) Il==inf ==> inf init slope = 0
fun.logpolynomial3.fun=@(p,I) p(2).*log(1+I).^p(1)+p(3);
fun.logpolynomial3.nrpara=3;
fun.logpolynomial3.LB=[0 0 0];
fun.logpolynomial3.X0=[1 1 0];
fun.logpolynomial3.UB=[10 10 10];
fun.logpolynomial3.var=2;

% Il=0 => p(3) Il=inf  =>0 bend at 10^p(1) slope at bend =  -p(3)/(4*p(2))/dec
fun.symLogisticsloginv.fun=@(p,I) (p(3)./(exp(-p(1)/p(2)).*I.^(1./(p(2)*log(10)))+1));
fun.symLogisticsloginv.nrpara=3;
fun.symLogisticsloginv.LB=[-10,0,0];
fun.symLogisticsloginv.X0=[1,1,0.5];
fun.symLogisticsloginv.UB=[10,20,1];
fun.symLogisticsloginv.var=2;

% Il=0 => 1 Il=inf  =>0 bend at 10^p(3)and 10^p(5) slope at bend1 =
% -p(2)/(4*p(4))/dec at bend1 =  -(1-p(2))/(4*p(6))/dec
fun.doublesymLogisticslogmina.fun =...
    @(p,I) p(1).*(1-p(2)./(1+exp(p(3)./p(4)).*I.^(-1./(p(4)*log(10))))-(1-p(2))./(1+exp(p(5)./p(6)).*I.^(-1./(p(6)*log(10)))));
    %@(p,Il)  p(1).*(1-p(2)./(1+exp(-(log10(Il)-p(3))/(p(4))))-(1-p(2))./(1+exp(-(log10(Il)-p(5))/p(6))));
fun.doublesymLogisticslogmina.nrpara = 6;
fun.doublesymLogisticslogmina.LB = [0,0,-10,0,-10,0];
fun.doublesymLogisticslogmina.X0 = [1,0.5,0,1/8,3,1/2];%  first bend at Il = 1 & -1/dec second bend at  Il = 1000 -1/4/dec
fun.doublesymLogisticslogmina.UB = [10,1,10,20,10,20];
fun.doublesymLogisticslogmina.var = 2; % function for Il or V or both

% Il=0 => 1 Il=inf  =>0 bend at 10^p(3)and 10^p(5) slope at bend1 =
% -p(2)/(4*p(4))/dec at bend1 =  -(1-p(2))/(4*p(6))/dec
fun.doublesymLogisticslogmina7.fun =...
    @(p,I) p(1).*(1-p(2)./(1+exp(p(3)./p(4)).*I.^(-1./(p(4)*log(10))))-(1-p(2))./(1+exp(p(5)./p(6)).*I.^(-1./(p(6)*log(10)))))+p(7);
    %@(p,Il)  p(1).*(1-p(2)./(1+exp(-(log10(Il)-p(3))/(p(4))))-(1-p(2))./(1+exp(-(log10(Il)-p(5))/p(6))));
fun.doublesymLogisticslogmina7.nrpara = 7;
fun.doublesymLogisticslogmina7.LB = [0,0,-10,0,-10,0,0];
fun.doublesymLogisticslogmina7.X0 = [1,0.5,0,1/8,3,1/2,0];%  first bend at Il = 1 & -1/dec second bend at  Il = 1000 -1/4/dec
fun.doublesymLogisticslogmina7.UB = [10,1,10,20,10,20,1];
fun.doublesymLogisticslogmina7.var = 2; % function for Il or V or both

% Il=0 => p(3)+p(4) Il=inf  =>0 bend at 10^p(1) slope at bend =  -p(3)/(4*p(2))/dec
fun.symLogisticsloginvd0.fun=@(p,I) (p(3)./(exp(-p(1)/p(2)).*I.^(1./(p(2)*log(10)))+1))+double(I==0)*p(4);
fun.symLogisticsloginvd0.nrpara=4;
fun.symLogisticsloginvd0.LB=[-10,0,0,0];
fun.symLogisticsloginvd0.X0=[1,1,0.5,0];
fun.symLogisticsloginvd0.UB=[10,20,1,10];
fun.symLogisticsloginvd0.var=2;


fun.weibull3para.fun=@(p,Il) p(3).*exp(-(Il./p(1)).^p(2));
fun.weibull3para.nrpara=3;
fun.weibull3para.LB=[0,0,0];
fun.weibull3para.X0=[50,0.1,0.4];
fun.weibull3para.UB=[1e2,10,1e2];
fun.weibull3para.var=2;

fun.weibull.fun=@(p,I) 1-exp(-(I./p(1)).^p(2));
fun.weibull.nrpara=2;
fun.weibull.LB=[0 0];
fun.weibull.X0=[1 1];
fun.weibull.UB=[1e5,1e5];
fun.weibull.var=2;

% Il=0 ==> p(1) Il = inf ==> 0 bend depends on value slope is asym
fun.asymLogisticsloginv.fun=@(p,Il) p(1).*(1-1./(1+exp(p(2)./p(3)).*Il.^(-1./(p(3)*log(10)))+exp(p(4)./p(5)).*Il.^(-1./(p(5)*log(10)))));
fun.asymLogisticsloginv.nrpara=5;
fun.asymLogisticsloginv.LB=[0,-10,0,-10,0];
fun.asymLogisticsloginv.X0=[1,0,1,3,1];
fun.asymLogisticsloginv.UB=[10,10,20,10,20];
fun.asymLogisticsloginv.var=2;

% Il=0 => 0 Il=inf  => p(1) depend on p(3) bend at p(2) slope at bend =  p(1)/(4*p(2))
fun.Logistics.fun=@(p,V) p(1)./(1+exp(-(V-p(2))./p(3)));
fun.Logistics.nrpara = 3;
fun.Logistics.LB=[0 -100 -1000];
fun.Logistics.X0=[1 -50 10];
fun.Logistics.UB=[100 100 1000];
fun.Logistics.var=1;

% Il=0 => p(4) Il=inf  => p(1) depend on p(3) bend at p(2) slope at bend =  p(1)/(4*p(2))
fun.Logistics4.fun=@(p,V) p(1)./(1+exp(-(V-p(2))./p(3)))+p(4);
fun.Logistics4.nrpara=4;
fun.Logistics4.LB=[0 -100 -1000,0];
fun.Logistics4.X0=[1 -50 10,0];
fun.Logistics4.UB=[100,100,1000,10];
fun.Logistics4.var=1;

% sum of two products of logistics for Il and V
fun.TauI.fun=@(p,X) p(1)./(1+p(2).*X(:,2).^p(3)).*(p(4)./(1+exp(-(X(:,1)-p(5))./p(6)))+p(7))+p(8)./(1+p(9).*X(:,2).^p(10)).*(p(11)./(1+exp((X(:,1)-p(12))./p(13)))+p(14));
fun.TauI.nrpara=14;
fun.TauI.LB=zeros(1,14);fun.TauI.LB(1,[5,12])=[-100,-100];
fun.TauI.X0=[1.738,1000,3,1,-20,20,1,0.26,0.01,2,1,-20,20,1];
fun.TauI.UB=[10,1e5,10,100,100,100,100,1,1,10,100,100,100,100];
fun.TauI.var=[1,2];

%0? not sure if possible
fun.One.fun=@(p,V) 1;
fun.One.nrpara=0;
fun.One.LB=[];
fun.One.X0=[];
fun.One.UB=[];
fun.One.var=1;

% rectification of conductance % be awere to correct p(1) depending on
% units used in target
fun.Grect.fun=@(p,V) p(1).*(1-p(2).*exp(-V./p(3)));
fun.Grect.nrpara=3;
fun.Grect.LB=[0,1.1,0];
fun.Grect.X0=[10,10,50];
fun.Grect.UB=[1e2,1e2,5*1e2];
fun.Grect.var=1;

fun.Grect_cont.fun=@(p,V) p(1).*(1-exp(-(V-p(2))./p(3)));
fun.Grect_cont.nrpara=3;
fun.Grect_cont.LB=[0,-100,0];
fun.Grect_cont.X0=[10,0,50];
fun.Grect_cont.UB=[1e2,100,5*1e2];
fun.Grect_cont.var=1;

% normal HH type
fun.GHH.fun=@(p,V) p(1).*(V-p(2));
fun.GHH.nrpara=2;
fun.GHH.LB=[0 -100];
fun.GHH.X0=[30 0];
fun.GHH.UB=[100 100];
fun.GHH.var=1;


%test functions based on alpha and beta

fun.OinfIV.fun = @(p,X) 1./(1+(1+exp(-(X(:,1)-p(3))./p(4)))./p(5).*exp(p(1)./p(2)).*X(:,2).^(-1./(p(2)*log(10))));
fun.OinfIV.nrpara = 5;
fun.OinfIV.LB = [-10 0 -100 -1000 0];
fun.OinfIV.X0 = [1 1 -50 10 1];
fun.OinfIV.UB = [10 20 100 1000 10];
fun.OinfIV.var = [1,2];

fun.TauOIV.fun = @(p,X) p(5)./(1+exp(-(X(:,1)-p(3))./p(4))+p(5).*exp(-p(1)./p(2)).*X(:,2).^(1./(p(2)*log(10))));
fun.TauOIV.nrpara = 5;
fun.TauOIV.LB = [-10 0 -100 -1000 0]; %boundaries are based on Target in s
fun.TauOIV.X0 = [1 1 -50 10 0.5];
fun.TauOIV.UB = [10 20 100 1000 10];
fun.TauOIV.var = [1,2];

fun.TauDAIV.fun = @(p,X) p(3).*p(7)./(p(3).*p(7).*exp(-p(1)./p(2)).*X(:,2).^(1./(p(2)*log(10)))+(1+exp(-(X(:,1)-p(5))./p(6))).*...
    (exp(-p(1)./p(2)).*(1-p(3)).*X(:,2).^(1./(p(2)*log(10)))+exp((p(4)-p(1))./(p(2)))));
fun.TauDAIV.nrpara = 7;
%                 p1  p2  p3    p4   p5     p6   p7
fun.TauDAIV.LB = [-10, 0,  0,   -10, 100, -1000, 0]; %boundaries are based on Target in s
fun.TauDAIV.X0 = [1,   1,  0.5, 1,   10,  10,    1];
fun.TauDAIV.UB = [10,  20, 1,   10,  100, 1000,  100];
fun.TauDAIV.var = [1,2];

fun.DAinfIV.fun = @(p,X) 1-(p(2).*p(6))./(p(2).*p(6)+(1-p(2)).*(1+exp(-(X(:,1)-p(4))./p(5)))+...
    exp(p(3)./p(1)).*X(:,2).^(-1./(p(1)*log(10))));
fun.DAinfIV.nrpara = 6;
fun.DAinfIV.LB = [0,  0,   -10, -100, -1000, 0];
fun.DAinfIV.X0 = [1,  0.5, 1,   -50,  10,    1];
fun.DAinfIV.UB = [20, 1,   10,  100,  1000,  100];
fun.DAinfIV.var = [1,2];

% double Il below
fun.TauDAdIlIV.fun = @(p,X) p(8)./(p(8).*exp(-p(1)./p(2)).*X(:,2).^(1./(p(2)*log(10)))+(1+exp(-(X(:,1)-p(6))./p(7))).*...
    (exp(-p(1)./p(2)).*X(:,2).^(1./(p(2)*log(10)))+exp(-p(3)./p(4)).*X(:,2).^(1./(p(4)*log(10)))+p(5)));
fun.TauDAdIlIV.nrpara = 8;
%                 p1  p2  p3    p4   p5     p6   p7     p8
fun.TauDAdIlIV.LB = [-10, 0,  -10, 0,  0,   100, -1000, 0]; %boundaries are based on Target in s
fun.TauDAdIlIV.X0 = [1,   1,  1,   1,  1,   10,  10,    1];
fun.TauDAdIlIV.UB = [10,  20, 10,  20, 100, 100, 1000,  100];
fun.TauDAdIlIV.var = [1,2];

fun.DAinfdIlIV.fun = @(p,X) 1-(p(8)./(p(8)+(1+exp(-(X(:,1)-p(6))./p(7))).*(1+...
    exp(-p(3)./p(4)).*exp(p(1)./p(2)).*X(:,2).^((p(2)-p(4))./(p(2).*p(4)*log(10)))+p(5).*exp(p(1)./p(2)).*X(:,2).^(-1./(p(2)*log(10))))));
fun.DAinfdIlIV.nrpara = 8;
fun.DAinfdIlIV.LB = [-10, 0,  -10, 0,  0,   100, -1000, 0];
fun.DAinfdIlIV.X0 = [1,   1,  1,   1,  1,   10,  10,    1];
fun.DAinfdIlIV.UB = [10,  20, 10,  20, 100, 100, 1000,  100];
fun.DAinfdIlIV.var = [1,2];

if newfunctions_flag
    fun = adjFields(fun,newfun,2);
end



%select names
% combination functions
prodfun = @(funV,funI) funV.*funI;
recsumfun = @(funV,funI) ((funV).^(-1)+(funI).^(-1)).^(-1);

% sets with multiple possibilities will probably not work anymore. Current
% version requires names cells to have same length. However! some
% optimization is required as via this method it is possible that during
% the intermediate fit, the same dependencies are fit multiple times.!
switch set
    case 'All1'
        OinfIfun_names = {'symLogisticslog','logarithm'};
        TauOIfun_names = {'symLogisticsloginv'};
        DAinfIfun_names = {'symLogisticsloga','polynomial'};
        TauDAIfun_names = {'asymLogisticsloginv'};
        TauOVfun_names = {'Logistics'};
        TauDAVfun_names = {'Logistics'};
        Grectfun_names = {'Grect'};
        combTauO = prodfun;
        combTauDA = prodfun;
    case 'All2'
        OinfIfun_names = {'symLogisticslog','logarithm'};
        TauOIfun_names = {'symLogisticsloginv'};
        DAinfIfun_names = {'symLogisticsloga','polynomial'};
        TauDAIfun_names = {'asymLogisticsloginv'};
        TauOVfun_names = {'Logistics'};
        TauDAVfun_names = {'Logistics'};
        Grectfun_names = {'Grect'};
        combTauO = prodfun;
        combTauDA = prodfun;
    case 'All3'
        OinfIfun_names = {'symLogisticslog','logarithm'};
        TauOIfun_names = {'symLogisticsloginv'};
        DAinfIfun_names = {'symLogisticsloga','polynomial'};
        TauDAIfun_names = {'TauI'};
        TauOVfun_names = {'Logistics'};
        TauDAVfun_names = {'None'};
        Grectfun_names = {'Grect'};
        combTauO = prodfun;
        combTauDA = prodfun;
    case 'standard'
        OinfIfun_names = {'symLogisticslog'};
        TauOIfun_names = {'symLogisticsloginv'};
        DAinfIfun_names = {'symLogisticslogmina'};
        TauDAIfun_names = {'doublesymLogisticslogmina'};
        TauOVfun_names = {'Logistics'};
        TauDAVfun_names = {'Logistics'};
        Grectfun_names = {'GHH'};
        combTauO = prodfun;
        combTauDA = prodfun;
    case 'IVfixed'
        OinfIfun_names = {'OinfIV'};
        TauOIfun_names = {'TauOIV'};
        DAinfIfun_names = {'DAinfIV'};
        TauDAIfun_names = {'TauDAIV'};
        TauOVfun_names = {'One'};
        TauDAVfun_names = {'One'};
        Grectfun_names = {'Grect'};
        combTauO = prodfun;
        combTauDA = prodfun;
    otherwise
        error('wrong input')
end

if newnames_flag
    if any(strcmpi(varargin,'OinfI'))
        OinfIfun_names = varargin{find(strcmpi(varargin,'OinfI'))+1};
    end
    if any(strcmpi(varargin,'TauOI'))
        TauOIfun_names = varargin{find(strcmpi(varargin,'TauOI'))+1};
    end
    if any(strcmpi(varargin,'DAinfI'))
        DAinfIfun_names = varargin{find(strcmpi(varargin,'DAinfI'))+1};
    end
    if any(strcmpi(varargin,'TauDAI'))
        TauDAIfun_names = varargin{find(strcmpi(varargin,'TauDAI'))+1};
    end
    if any(strcmpi(varargin,'TauOV'))
        TauOVfun_names = varargin{find(strcmpi(varargin,'TauOV'))+1};
    end
    if any(strcmpi(varargin,'TauDAV'))
        TauDAVfun_names = varargin{find(strcmpi(varargin,'TauDAV'))+1};
    end
    if any(strcmpi(varargin,'Grect'))
        Grectfun_names = varargin{find(strcmpi(varargin,'Grect'))+1};
    end
    if any(strcmpi(varargin,'combTauO'))
        combTauO = varargin{find(strcmpi(varargin,'combTauO'))+1};
    end
    if any(strcmpi(varargin,'combTauDA'))
        combTauDA = varargin{find(strcmpi(varargin,'combTauDA'))+1};
    end
end
end