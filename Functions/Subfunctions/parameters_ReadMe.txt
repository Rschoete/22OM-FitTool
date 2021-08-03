parameters in fit22HH

HPC_flag			if HPC_flag: all display and plots off!
time_flag: 			record time or not
display_flag: 		display info
newfun_flag:		new function to add to databaqe
plotFeat_flag:		flag plot features
plotTau_flag:		flag plot Tau
plot_flag
Colors: 			function that allows for iteration over specific set of collors
type:				function allows iteration over different linetypes

checkDouble_flag: 	check for doubles in target files and combine them
featNames: 			feature names for which checkDouble applies (in future version mayby also indicative for featurenames to fit on)
featNames_SD: 		same as above but now standard deviation values
featNamesTr: 		featnames tau recovery has to be changed!!
featNamesTr_SD:		same above also

sort_flag: 			sort values in target or not
plotFeat_flag: 		create figures of the features and their depence to Il and V
featNamesT_plot:	names of Taufeatures to plot
featNamesI_plot:	names of I features to plot

errorbar_flag:		plot with errorbar or not
fig_pos: 			positions where figures are to be plotted

maxIter:			preset maximum iterations
maxFunEval:			preset maximum function iterations
DisplayIter: 		display results during iterations,
					'final','iter','off'
FunTol				Tolerance on function values for considering solutions equal, specified as a nonnegative scalar.
					Solvers consider two solutions identical if they are within XTolerance relative distance of each other and
					have objective function values within FunctionTolerance relative difference of each other
XTol				X tolerance in multistart

IVdepend_targets:	first fit based on 'features' or 'postcurrenttraces'
	(Oinf, TauO, DAinf, tauDA,g) already extracted and in Target file, 'currenttraces'
	exract aforementioned here
recal_opt_FBIVdepend
		.TauO
			.dt		Input not necessary
			.tend	Input not necessary
			.optimoptions
				{'lsqcurvefit','Algorithm','trust-region-reflective','MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter}
			.lb		Input not necessary
			.ub		Input not necessary
			.msopt multistart options
				{'Display',DisplayIter,'UseParallel','always','FunctionTolerance',FunTol,'XTolerance',XTol};
			.msnr multistart nrstart
		.TauDA
			.dt		Input not necessary
			.tend	Input not necessary
			.ub		Input not necessary
			.ubr (upper boundary of recovery param (higher than normal ub)) Input not necessary
			.optimoptions
				{'lsqcurvefit','Algorithm','trust-region-reflective','MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter}
			.lb	Input not necessary
			.msopt
				{'Display',DisplayIter,'UseParallel','always','FunctionTolerance',funTol,'XTolerance',funTol};
			.msnr
OinfDAinf_flag: 	if above method featurebased Oinf and DAinf need to be extracted still,
		based on Ipeak, Iss ==> value = 0 (no Oinf and DA yet)

%Load functionDB
set: 		which of predefined sets should be used? inputs: All1, All2, All3, standard (default)
newfun: 	structure of with new functions or old functions but changed parameters (eg LB,UB,XO)
		==> if included activates newfunctions_flag
newnames_flag: which functions are tested activated when
		varargin contains 'OinfI', 'TauOI', 'DAinfI','TauDAI','TauOV','TauDAV','Grect',combTauO,combTauDa


%fitTauO
options_TauO
	.optimoptions
		 {'lsqcurvefit','Algorithm','trust-region-reflective','MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
	.msopt
		{'Display',DisplayIter,'UseParallel','always','FunctionTolerance',funTol,'XTolerance',funTol};
	.msnr

%fitTauDA
options_TauDA
	.optimoptions
		{'lsqcurvefit','Algorithm','trust-region-reflective','MaxFunctionEvaluations',600e3,'MaxIterations',maxIter};
	.msnr
	.msopt
		{'Display',Displayiter,'UseParallel','always','FunctionTolerance',1e-10,'XTolerance',1e-10};
%fit Oinf DAinf, Grect
wGODA.type
wGoda.w  = 1			size = [target_sp,2] left weihgts Ipeak, right weights Iss
options_GODA
	.I_nonlcon			I values to included in nonlinconstraint
	.V1_nonlcon 		Values for which Grect has to be positive
	.V2_nonlcon			Values for which Grect has to be negative
	.extrac_nonlcon		Extraconstraints
		.IO_smaller 	Intensities extra conditions Oinf
		.Olim_smaller 	limits of which O has to be smaller
		.IO_bigger 		Intensities extra conditions Oinf
		.Olim_bigger 	limits of which O has to be smaller
		.IDA_smaller 	Intensities extra conditions DAinf
		.DAlim_smaller 	limits of which DA has to be smaller
		.IDA_bigger 	Intensities extra conditions DAinf
		.DAlim_bigger 	limits of which DA has to be smaller
	.extraceq_nonlcon
		.VG 			value when Grect ==0

	.optimoptions
		{'fmincon','Algorithm','sqp','MaxFunctionEvaluations',maxFunEval,'MaxIterations',maxIter};
	.msopt
		{'Display',DisplayIter,'UseParallel','always','FunctionTolerance',funTol,'XTolerance',XTol};
	.msnr
	.dt
	.tend

% run ODE
	OdeOpts
		{'RelTol',1e-6,'AbsTol',1e-6,'Maxstep',0.001};

% final optimization

fopt_select:		'all' or 'best'  final fit again over all possible combinations or only the best
fopt_objMethod:		'currenttraces' or 'features' cost function is based on features or currenttraces
fopt_solverMethod:	'ms', 'pso','ga' only ms encoded for now
options_fopt
	.incl_tr_flag   include tau recov or not
	.fLB			factor for lowerbounds (multiplied with X0)
	.fUB			factor for upperbounds

	.CT				only necessary if objMethod=currenttraces
		.intervals
		.method_CT	this is about scaling obj function
					'normalized' (with peak)
					'zscore'
					'weighted'
					'normalized+weighted'
		.w			needed when above weighted assigned
					has to be a column or matrix with number of rows equal to number of targets tV (take care this can change when doublecheck is on)
					also sort flag affects this ==> encode autofix!
		.wtr		same as above but now for tau recovery
		.obj		'rms' or 'q2'

	.FB				only necessary when objMethod features
		.tend_sp	end time whereoff features extracted single pulse necessary if no time field in target
		.tend_tr
		.method_tr	how extract Taurecov feature
					'estimpowerb'
					'noestimpowerb'
					'currenttraces'
		.trestim 	if estimpowerb selected
			.dt		optional
			.optimoptions
			.lb		optionall
			.ub		optionall
			.msopt
			.msnr
		.intervals

		.dt_tr		needed when currenttraces
		.xx  		estim for interpolation to TauRecov
		.dt
		.nrr_featextr	nrruns in multistart
		.plot_featextr
		.varargin{
		X0			has to be column or matrix with length tV
		LB
		UB
		}
		.w			matrix 6columns and tV rows

	.I_nonlcon			I values to included in nonlinconstraint
	.V1_nonlcon 		Values for which Grect has to be positive
	.V2_nonlcon			Values for which Grect has to be negative
	.extrac_nonlcon		Extraconstraints
		.IO_smaller 	Intensities extra conditions Oinf
		.Olim_smaller 	limits of which O has to be smaller
		.IO_bigger 		Intensities extra conditions Oinf
		.Olim_bigger 	limits of which O has to be smaller
		.IDA_smaller 	Intensities extra conditions DAinf
		.DAlim_smaller 	limits of which DA has to be smaller
		.IDA_bigger 	Intensities extra conditions DAinf
		.DAlim_bigger 	limits of which DA has to be smaller
	.extraceq_nonlcon
		.VG 			value when Grect ==0
	.Outfun
		.OdeOpts
	.ms
		.optimoptions
		.msopt
		.msnr