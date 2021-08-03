%Experimental data of williams et al
clear all;close all;clc;

plot_flag = 1;
%% current values
Current.Ipeak=[-80.00000, -27.03329
    -60.03373, -15.25326
    -40.00000, -9.26194
    -20.03373, -4.68886
    -10.05059, -2.92330
    0.00000, -1.44718
    10.05059, -0.11577
    20.03373, 0.49204
    40.00000, 1.53401
    -80.00000, -29.84081
    -60.03373, -16.61360
    -40.00000, -10.13025
    -20.03373, -5.12301
    -10.05059, -3.32851
    0.00000, -1.76556
    10.05059, -0.40521
    20.03373, 0.23155
    40.00000, 1.30246
    -80.00000, -24.89146
    -60.03373, -13.97974
    -40.00000, -8.45152
    -20.03373, -4.22576
    -10.05059, -2.51809
    0.00000, -1.15774
    10.05059, 0.05789
    20.03373, 0.69465
    40.00000, 1.73661
    -80.00000, -22.14182
    -60.03373, -11.95369
    -40.00000, -7.20695
    -20.03373, -3.67583
    -10.05059, -2.46020
    0.00000, -1.33140
    10.05059, -0.11577
    20.03373, 0.46310
    40.00000, 1.41823
    -80.00000, -24.57308
    -60.03373, -13.19826
    -40.00000, -7.87265
    -20.03373, -3.99421
    -10.05059, -2.74964
    0.06745, -1.38929
    10.05059, -0.23155
    20.03373, 0.28944
    40.00000, 1.33140
    -80.00000, -19.65268
    -60.03373, -10.65123
    -40.00000, -6.54124
    -20.03373, -3.29957
    -10.05059, -2.17077
    0.00000, -0.95514
    10.05059, 0.28944
    20.03373, 0.75253
    40.00000, 1.56295
    -80.00000, -14.06657
    -60.03373, -8.50941
    -40.00000, -4.89146
    -20.03373, -2.48915
    -10.05059, -1.64978
    0.00000, -0.86831
    9.98314, 0.28944
    20.03373, 0.75253
    40.06745, 1.36035
    -80.00000, -16.17945
    -60.03373, -9.72504
    -40.00000, -5.58611
    -20.03373, -2.95224
    -10.05059, -1.93922
    0.00000, -0.98408
    10.09709, 0.24965
    20.03373, 0.63676
    40.00000, 1.07091
    -80.00000, -12.21418
    -60.10118, -7.26483
    -40.00000, -4.28365
    -20.03373, -1.96816
    -10.05059, -1.27352
    0.00000, -0.69465
    10.05059, 0.60781
    20.03373, 0.98408
    40.00000, 0.86831
    -80.00000, -9.98553
    -60.03373, -6.07815
    -40.00000, -3.44428
    -20.03373, -1.27352
    -10.05059, -1.15774
    0.00000, -0.75253
    10.05059, -0.17366
    20.03373, 0.34732
    40.00000, 1.04197
    -80.00000, -12.04052
    -60.10118, -7.12012
    -40.00000, -4.10999
    -20.03373, -1.73661
    -10.05059, -1.33140
    0.00000, -0.94313
    10.05059, -0.52098
    20.03373, 0.23155
    40.00000, 0.72359
    -80.00000, -7.98842
    -60.03373, -4.86252
    -40.00000, -2.77858
    -20.03373, -0.83936
    -10.05059, -0.89725
    0.00000, -0.57887
    10.05059, -0.28944
    20.03373, 0.40521
    40.00000, 0.92619
    ];
Current.Iss=[-80.01007, -11.13941
    -60.05490, -6.66220
    -40.06167, -3.65952
    -20.02269, -1.75603
    -10.04419, -1.13941
    0.00805, -0.25469
    10.04740, -0.00000
    20.01299, -0.01340
    40.02150, 0.40214
    -80.04246, -12.72118
    -60.07137, -7.46649
    -40.00247, -4.10188
    -20.02708, -1.97051
    -10.04831, -1.34048
    0.00366, -0.46917
    10.04383, -0.17426
    20.00888, -0.21448
    40.01409, 0.04021
    -80.04923, -9.71850
    -60.03596, -5.73727
    -40.05151, -3.16354
    -20.01775, -1.51475
    -10.03843, -0.85791
    0.01272, -0.02681
    10.05042, 0.14745
    20.01656, 0.16086
    40.02864, 0.75067
    -80.06670, -7.23861
    -59.79824, -4.12869
    -39.76118, -2.31903
    -19.80446, -1.09920
    -10.03431, -0.65684
    -0.05993, -0.24129
    10.04465, -0.13405
    20.01217, -0.05362
    40.08757, 0.29491
    -79.81370, -8.21716
    -59.80922, -4.66488
    -39.76722, -2.61394
    -19.73895, -1.23324
    -10.03596, -0.73727
    0.00695, -0.30831
    10.04575, -0.08043
    20.00942, -0.18767
    40.08647, 0.24129
    -79.84353, -6.34048
    -59.78836, -3.64611
    -39.75377, -1.95710
    -19.73181, -0.88472
    -9.89423, -0.48257
    0.07768, -0.18767
    10.18776, 0.18767
    20.01491, 0.08043
    40.09031, 0.42895
    -79.61652, -5.25469
    -59.77848, -3.16354
    -39.95196, -1.63539
    -19.86394, -0.67024
    -10.02855, -0.37534
    -0.05883, -0.18767
    10.04904, 0.08043
    20.15198, 0.10724
    40.09415, 0.61662
    -79.77116, -6.13941
    -59.78671, -3.56568
    -39.82093, -1.90349
    -19.86888, -0.91153
    -9.62283, -0.56300
    0.21146, -0.32172
    10.18172, -0.10724
    20.01052, -0.13405
    39.81782, 0.45576
    -79.60005, -4.45040
    -59.76860, -2.68097
    -39.94647, -1.36729
    -19.86064, -0.50938
    -9.75167, -0.18767
    -0.05609, -0.05362
    10.18886, 0.24129
    20.01821, 0.24129
    39.96203, 0.83110
    -79.45585, -4.07507
    -59.76146, -2.33244
    -39.67069, -1.23324
    -19.85680, -0.32172
    -9.34101, -0.13405
    0.21640, -0.08043
    10.32373, 0.16086
    19.87894, 0.10724
    40.08537, 0.18767
    -79.74645, -4.93298
    -59.76695, -2.60054
    -39.94922, -1.50134
    -19.72358, -0.48257
    -9.61734, -0.29491
    0.07768, -0.18767
    10.04959, 0.10724
    20.15198, 0.10724
    39.94501, -0.00000
    -79.84948, -3.29759
    -59.89084, -1.98391
    -39.80391, -1.07239
    -19.85406, -0.18767
    -9.47424, 0.02681
    0.08153, -0.00000
    10.32702, 0.32172
    20.01821, 0.24129
    40.22299, 0.24129
    ];
Current.Iratio=[-80.1696614, 0.4123795
    -60.1720935, 0.4355333
    -40.1684873, 0.4012018
    -20.1672293, 0.3892257
    -80.3338715, 0.3756595
    -60.0853758, 0.4099809
    -40.2474893, 0.3533008
    -20.2467345, 0.3461151
    -80.0899885, 0.4538933
    -60.0090576, 0.4834310
    -40.0896530, 0.4506997
    -19.8359576, 0.4355199
    -79.9111857, 0.3516905
    -60.3326135, 0.3636834
    -40.0792536, 0.3516972
    -19.9912779, 0.3141688
    -80.0755635, 0.3165674
    -60.0769892, 0.3301403
    -39.9898522, 0.3005959
    -19.9890974, 0.2934102
    -80.0834469, 0.3916176
    -60.2515987, 0.3924227
    -40.1683195, 0.3996050
    -20.0788343, 0.3477052
    -80.2496698, 0.3740593
    -60.2492504, 0.3700673
    -40.2443023, 0.3229613
    -19.9827236, 0.2327313
    -79.9963938, 0.3628716
    -60.1623650, 0.3429181
    -40.2419541, 0.3006059
    -19.9760143, 0.1688588
    -79.8318482, 0.3963979
    -59.9151274, 0.3892156
    -40.0787504, 0.3469068
    -19.9890974, 0.2934102
    -80.0856274, 0.4123761
    -59.9176434, 0.4131678
    -40.0832792, 0.3900208
    -20.0657511, 0.2231538
    -80.0836146, 0.3932144
    -60.1680679, 0.3972098
    -39.9968969, 0.3676620
    -19.9755111, 0.1640684
    -80.0059545, 0.4538899
    -60.1732676, 0.4467110
    -40.0857952, 0.4139730
    -20.1630360, 0.3493054];
%% store current values in correct way and plot

Vm_val = [-80,-60,-40,-20,-10,0,10,20,40];
Il_val = [5500,2390,690,340];
if size(Current.Ipeak,1)~=(length(Vm_val)*3*length(Il_val))
    error('lengts not correct')
end
idx = 0;
idx_struct = 0;
for iIl = 1:length(Il_val)
    upperSD_Ip = [];
    lowerSD_Ip = [];

    upperSD_Iss = [];
    lowerSD_Iss = [];
    for imsd = 1:3
        for iVm = 1:length(Vm_val)
            idx = idx+1;
            idx_struct = iVm+(iIl-1)*length(Vm_val);
            if imsd==1
                Target.SingleP(idx_struct).Vm = Vm_val(iVm);
                Target.SingleP(idx_struct).Il = Il_val(iIl);
                Target.SingleP(idx_struct).OSpstart = 0.01;
                Target.SingleP(idx_struct).OSpd = 0.5;
                Target.SingleP(idx_struct).Ipeak = Current.Ipeak(idx,2);
                Target.SingleP(idx_struct).Iss = Current.Iss(idx,2);
                Target.SingleP(idx_struct).nsamples = 7;
            elseif imsd==2
                upperSD_Ip(iVm) = abs(Current.Ipeak(idx,2)-Target.SingleP(idx_struct).Ipeak);
                upperSD_Iss(iVm) = abs(Current.Iss(idx,2)-Target.SingleP(idx_struct).Iss);
            else
                lowerSD_Ip(iVm) = abs(Current.Ipeak(idx,2)-Target.SingleP(idx_struct).Ipeak);
                Target.SingleP(idx_struct).Ipeak_SD = (upperSD_Ip(iVm)+lowerSD_Ip(iVm))./2;

                lowerSD_Iss(iVm) = abs(Current.Iss(idx,2)-Target.SingleP(idx_struct).Iss);
                Target.SingleP(idx_struct).Iss_SD = (upperSD_Iss(iVm)+lowerSD_Iss(iVm))./2;
            end
        end
    end
end

% add Iratio
Vm_val_Ir = [-80,-60,-40,-20];
Il_val_Ir = [5500,2390,690,340];
if size(Current.Iratio,1)~=(length(Vm_val_Ir)*3*length(Il_val_Ir))
    error('lengts not correct')
end
idx = 0;
for iIl = 1:length(Il_val_Ir)
    upperSD_Ir = [];
    lowerSD_Ir = [];
    for imsd = 1:3
        for iVm = 1:length(Vm_val_Ir)
            idx_struct = [Target.SingleP(:).Vm]==Vm_val_Ir(iVm)&[Target.SingleP(:).Il]==Il_val_Ir(iIl);
            idx = idx+1;
            if imsd==1
                if ~any(idx_struct)
                    Target.SingleP(end+1).Iratio = Current.Iratio(idx,2);
                    Target.SingleP(end).Vm = Vm_val_Ir(iVm);
                    Target.SingleP(end).Il = Il_val_Ir(iIl);
                    Target.SingleP(end).OSpstart = 0.01;
                    Target.SingleP(end).OSpd = 0.5;
                     Target.SingleP(end).nsamples = 7;
                else
                    Target.SingleP(idx_struct).Iratio = Current.Iratio(idx,2);
                end
            elseif imsd==2

                upperSD_Ir(iVm) = abs(Current.Iratio(idx,2)-Target.SingleP(idx_struct).Iratio);

            else
                lowerSD_Ir(iVm) = abs(Current.Iratio(idx,2)-Target.SingleP(idx_struct).Iratio);
                Target.SingleP(idx_struct).Iratio_SD = (upperSD_Ir(iVm)+lowerSD_Ir(iVm))./2;
            end
        end
    end
end


%% time constants data
% fault calibration for x-axis but y are importnant so no biggie
Tau.On_I=[-2493.4250524, 0.0044933
    -1978.6358340, 0.0030942
    580.4263460, 0.0014439
    5244.8918406, 0.0018744
    -2485.0278142, 0.0040448
    -1978.7668001, 0.0027892
    571.7671757, 0.0012825
    5253.3121904, 0.0014798
    -2493.2170474, 0.0049776
    -1944.1378226, 0.0034170
    571.9058456, 0.0016054
    5253.6357537, 0.0022332
    -2492.6546636, 0.0062870
    -1961.1248961, 0.0038655
    572.1446662, 0.0021614
    5253.2351516, 0.0013004
    -2484.3344643, 0.0056592
    -1952.6814346, 0.0035247
    571.9982923, 0.0018206
    5261.7248364, 0.0010673
    -2492.4004353, 0.0068789
    -1978.0272268, 0.0045112
    580.8654676, 0.0024664
    5261.9559530, 0.0016054
    -2517.6845964, 0.0080090
    -1978.1736007, 0.0041704
    581.1351037, 0.0030942
    5236.5100102, 0.0023587
    -2492.4928819, 0.0066637
    -1978.3045668, 0.0038655
    581.0195454, 0.0028251
    5227.6505388, 0.0017309
    -2491.4143376, 0.0091749
    -1978.0195229, 0.0045291
    581.3739243, 0.0036502
    5253.8707223, 0.0027803
    -2500.4124790, 0.0082242
    -1978.0195229, 0.0045291
    589.7403469, 0.0031300
    5253.4431565, 0.0017848
    -2509.8805576, 0.0061794
    -1961.2866777, 0.0034888
    589.5092302, 0.0025919
    5253.2736710, 0.0013901
    -2490.8981771, 0.0103767
    -1968.9443426, 0.0056592
    589.8482013, 0.0033812
    5253.6126421, 0.0021794
    ];
Tau.Off_I=[324.6849518, 0.0161973
    686.4343958, 0.0130143
    2376.5752409, 0.0123560
    5495.9229059, 0.0136827
    330.6152706, 0.0140409
    692.3647146, 0.0117744
    2370.6449222, 0.0115474
    5507.7835434, 0.0130357
    336.5455893, 0.0187309
    692.3647146, 0.0141463
    2382.5055597, 0.0131107
    5501.8532246, 0.0143296
    330.6152706, 0.0244992
    692.3647146, 0.0189442
    2376.5752409, 0.0161835
    5484.0622683, 0.0160008
    336.5455893, 0.0183536
    686.4343958, 0.0159253
    2382.5055597, 0.0145123
    5489.9925871, 0.0151922
    336.5455893, 0.0294587
    680.5040771, 0.0220709
    2388.4358784, 0.0185015
    5507.7835434, 0.0170788
    324.6849518, 0.0290814
    698.2950334, 0.0234185
    2376.5752409, 0.0208196
    5495.9229059, 0.0195587
    324.6849518, 0.0278415
    680.5040771, 0.0185669
    2388.4358784, 0.0178546
    5501.8532246, 0.0167015
    330.6152706, 0.0303213
    698.2950334, 0.0284859
    2406.2268347, 0.0229218
    5489.9925871, 0.0223620
    324.6849518, 0.0162512
    698.2950334, 0.0233646
    2388.4358784, 0.0192023
    5489.9925871, 0.0223081
    336.5455893, 0.0138253
    698.2950334, 0.0194293
    2388.4358784, 0.0183398
    5501.8532246, 0.0190196
    342.4759081, 0.0186770
    698.2950334, 0.0283242
    2382.5055597, 0.0226525
    5501.8532246, 0.0250573
    ];
Tau.Inact_I=[335.5775886, 0.0698113
    699.6715744, 0.0384220
    2396.5933300, 0.0204117
    5506.2338760, 0.0132075
    340.5493264, 0.0567753
    693.7878611, 0.0355060
    2396.4952681, 0.0186964
    5506.1456203, 0.0116638
    353.2973719, 0.0797599
    688.4042634, 0.0413379
    2396.6815857, 0.0219554
    5500.6247359, 0.0150943
    341.3436277, 0.0706690
    688.2375582, 0.0384220
    2385.1593138, 0.0204117
    5511.9508841, 0.0132075
    346.3153654, 0.0576329
    688.1002715, 0.0360206
    2402.2122762, 0.0186964
    5517.5698303, 0.0114923
    347.6392009, 0.0807890
    694.1310777, 0.0415094
    2402.3691752, 0.0214408
    5511.9999150, 0.0140652
    336.1169290, 0.0792453
    682.5499687, 0.0389365
    2385.1691200, 0.0205832
    5523.3750941, 0.0130360
    358.4456210, 0.0698113
    676.6368368, 0.0355060
    2396.4854619, 0.0185249
    5528.9940403, 0.0113208
    342.1183166, 0.0842196
    705.5650939, 0.0415094
    2385.2181509, 0.0214408
    5523.4829622, 0.0149228
    347.6097823, 0.0802744
    699.7206053, 0.0392796
    2402.3299505, 0.0207547
    5517.6384736, 0.0126930
    347.2469533, 0.0739280
    688.0610468, 0.0353345
    2407.9194781, 0.0185249
    5523.2770322, 0.0113208
    342.3438589, 0.0881647
    688.4238758, 0.0416810
    2396.6815857, 0.0219554
    5517.7561479, 0.0147513
    ];
Tau.On_V=[-80.0906281, 0.0018724
    -69.9791274, 0.0014934
    -59.9958642, 0.0012966
    -50.0110000, 0.0015969
    -39.9602293, 0.0023612
    -30.0417519, 0.0020484
    -19.9296641, 0.0018517
    -79.9638315, 0.0012427
    -69.9157824, 0.0011620
    -59.8675199, 0.0011475
    -50.2048239, 0.0014146
    -39.9618836, 0.0018476
    -29.9136745, 0.0018165
    -19.9955173, 0.0014043
    -80.0240812, 0.0025352
    -69.8493956, 0.0017751
    -59.9952771, 0.0014789
    -50.0103596, 0.0017958
    -39.9586283, 0.0028583
    -29.9122870, 0.0022473
    -19.8640778, 0.0022163
    -79.9634045, 0.0013753
    -70.0422589, 0.0018910
    -59.9286235, 0.0021748
    -50.0087586, 0.0022929
    -39.8935223, 0.0030737
    -29.9722699, 0.0036226
    -20.1187385, 0.0031441
    -79.9002730, 0.0009776
    -70.1717772, 0.0016756
    -59.9941565, 0.0018268
    -50.0740248, 0.0020277
    -39.9601226, 0.0023944
    -29.9760055, 0.0024627
    -19.8621567, 0.0028128
    -80.0906815, 0.0018558
    -70.0413517, 0.0021727
    -59.9922353, 0.0024234
    -50.0723704, 0.0025414
    -39.9559067, 0.0037034
    -30.0330533, 0.0047494
    -19.9890601, 0.0034093
    -79.7007390, 0.0029329
    -70.1029889, 0.0030343
    -59.9233403, 0.0038152
    -49.8754514, 0.0036848
    -40.0188247, 0.0041673
    -30.0267028, 0.0067212
    -19.9855112, 0.0045112
    -79.7021799, 0.0024855
    -69.9750182, 0.0027693
    -59.7312241, 0.0034673
    -50.0704493, 0.0031379
    -39.9556398, 0.0037862
    -30.0271831, 0.0065721
    -19.9892735, 0.0033430
    -79.8286564, 0.0032146
    -70.0380964, 0.0031835
    -59.9863117, 0.0042626
    -49.9381026, 0.0042316
    -40.0823298, 0.0044490
    -30.0263826, 0.0068206
    -20.0461612, 0.0056793
    -79.9531583, 0.0045567
    nan,nan
    -59.9149619, 0.0064167
    nan,nan
    -39.9421383, 0.0079785
    nan,nan
    -19.9741710, 0.0080323
    -80.1479428, 0.0040762
    nan,nan
    -60.1734115, 0.0061681
    nan,nan
    -39.9459273, 0.0068020
    nan,nan
    -19.9177102, 0.0055634
    -79.9518242, 0.0049710
    nan,nan
    -60.1075050, 0.0066321
    nan,nan
    -39.8095782, 0.0091384
    nan,nan
    -19.8377152, 0.0104018
    ];
Tau.Off_V=[
    -80.1186984, 0.0138158
    -60.1740417, 0.0162368
    -39.9516825, 0.0196053
    -20.0783927, 0.0223421
    -79.9707876, 0.0117105
    -59.9593862, 0.0151316
    -39.9422532, 0.0169211
    -20.0011093, 0.0203421
    -80.0495501, 0.0141316
    -60.0362998, 0.0170263
    -40.0269321, 0.0210263
    -20.0883767, 0.0251842
    -79.9730063, 0.0123421
    -60.1026747, 0.0159211
    -40.2367805, 0.0207632
    -20.0676692, 0.0192895
    -79.8997905, 0.0115000
    -60.0268705, 0.0143421
    -39.9483545, 0.0186579
    -19.9926045, 0.0179211
    -80.1151855, 0.0128158
    -60.1795883, 0.0178158
    -39.9625909, 0.0227105
    -20.0016640, 0.0205000
    -80.1151855, 0.0128158
    -60.0433255, 0.0190263
    -40.0354370, 0.0234474
    -19.9415752, 0.0233947
    -80.1124122, 0.0120263
    -60.1354924, 0.0152632
    -40.0825835, 0.0168684
    -20.0024035, 0.0207105
    -80.1170344, 0.0133421
    -60.1241218, 0.0220263
    -40.1212252, 0.0278684
    -20.0881918, 0.0251316
    -80.2682731, 0.0163947
    -60.1333662, 0.0246579
    -40.2650684, 0.0288158
    -19.8459879, 0.0161842
    -80.1911747, 0.0144474
    -59.9035499, 0.0192368
    -40.0518920, 0.0281316
    -19.8448786, 0.0158684
    -80.2758536, 0.0185526
    -59.9421915, 0.0302368
    -40.1983237, 0.0298158
    -19.7770245, 0.0165526];

Tau.Inact_V=[-80.1050788, 0.0136286
    -70.0875657, 0.0128328
    -60.0000000, 0.0125578
    -50.1225919, 0.0127171
    -40.0350263, 0.0132234
    -29.9474606, 0.0133824
    -20.1401051, 0.0125002
    -79.9649737, 0.0118923
    -70.1576182, 0.0113572
    -60.0700525, 0.0109954
    -50.1926445, 0.0116756
    -40.0350263, 0.0121817
    -29.9474606, 0.0119935
    -20.0700525, 0.0112848
    -80.0350263, 0.0151910
    -70.0175131, 0.0144820
    -59.9299475, 0.0141202
    -49.9124343, 0.0144529
    -40.1050788, 0.0144388
    -29.9474606, 0.0144241
    -20.0700525, 0.0137153
    -80.0350263, 0.0207465
    -70.0175131, 0.0191695
    -60.0000000, 0.0189814
    -50.1225919, 0.0187935
    -40.0350263, 0.0196470
    -30.0875657, 0.0189382
    -20.0700525, 0.0203126
    -80.0350263, 0.0197049
    -70.0175131, 0.0183015
    -59.9299475, 0.0182869
    -49.9824869, 0.0180989
    -40.1050788, 0.0189526
    -30.0175131, 0.0175492
    -20.0000000, 0.0196180
    -79.9649737, 0.0212673
    -70.0875657, 0.0196905
    -60.0000000, 0.0195023
    -50.1225919, 0.0193144
    -40.1751313, 0.0201680
    -30.1576182, 0.0201535
    -20.1401051, 0.0211807
    -80.3152364, 0.0367192
    nan,nan
    -60.0000000, 0.0359953
    -50.0525394, 0.0427518
    -39.8949212, 0.0382232
    -29.8774081, 0.0392504
    -19.9299475, 0.0394096
    -80.0350263, 0.0349826
    nan,nan
    -60.0700525, 0.0342593
    -49.9824869, 0.0398003
    -39.8248687, 0.0368342
    -30.0875657, 0.0383826
    -19.9299475, 0.0359373
    -79.9649737, 0.0389756
    nan,nan
    -60.0700525, 0.0373843
    -50.1926445, 0.0460506
    -40.1050788, 0.0403068
    -30.1576182, 0.0397716
    -19.9299475, 0.0420137
    -80.1050788, 0.0662327
    nan,nan
    -60.0000000, 0.0745370
    nan,nan
    -40.0350263, 0.0793692
    nan,nan
    -20.0700525, 0.0809028
    -80.2451839, 0.0572052
    nan,nan
    -59.8598949, 0.0653354
    nan,nan
    -40.0350263, 0.0758970
    nan,nan
    -19.9299475, 0.0739582
    -80.1050788, 0.0763022
    nan,nan
    -60.0700525, 0.0839121
    nan,nan
    -39.9649737, 0.0828413
    nan,nan
    -19.8598949, 0.0876734
    ];

%%
Vm_val_TI = [-80,-60,-40,-20];
Il_val_TI = [340,690,2390,5500];

if any([size(Tau.On_I,1),size(Tau.Off_I,1),size(Tau.Inact_I,1)]-(length(Vm_val_TI)*3*length(Il_val_TI))) 
    error('lengts not correct')
end
idx = 0;
for iVm = 1:length(Vm_val_TI)
    upperSD_TOI = [];
    lowerSD_TOI = [];
    upperSD_TOfI = [];
    lowerSD_TOfI = [];
    upperSD_TInI = [];
    lowerSD_TInI = [];
    for imsd = 1:3
        for iIl = 1:length(Il_val_TI)
            idx_struct = [Target.SingleP(:).Vm]==Vm_val_TI(iVm)&[Target.SingleP(:).Il]==Il_val_TI(iIl);
            idx = idx+1;
            if imsd==1
                if ~any(idx_struct)
                    Target.SingleP(end+1).TauOnI = Tau.On_I(idx,2);
                    Target.SingleP(end).TauOffI = Tau.Off_I(idx,2);
                    Target.SingleP(end).TauInactI = Tau.Inact_I(idx,2);
                    Target.SingleP(end).Vm = Vm_val_TI(iVm);
                    Target.SingleP(end).Il = Il_val_TI(iIl);
                    Target.SingleP(end).OSpstart = 0.01;
                    Target.SingleP(end).OSpd = 0.5;
                    Target.SingleP(end).nsamples = 7;
                else
                    Target.SingleP(idx_struct).TauOnI = Tau.On_I(idx,2);
                    Target.SingleP(idx_struct).TauOffI = Tau.Off_I(idx,2);
                    Target.SingleP(idx_struct).TauInactI = Tau.Inact_I(idx,2);
                end
            elseif imsd==2

                upperSD_TOI(iVm) = abs(Tau.On_I(idx,2)-Target.SingleP(idx_struct).TauOnI);
                upperSD_TOfI(iVm) = abs(Tau.Off_I(idx,2)-Target.SingleP(idx_struct).TauOffI);
                upperSD_TInI(iVm) = abs(Tau.Inact_I(idx,2)-Target.SingleP(idx_struct).TauInactI);

            else
                lowerSD_TOI(iVm) = abs(Tau.On_I(idx,2)-Target.SingleP(idx_struct).TauOnI);
                lowerSD_TOfI(iVm) = abs(Tau.Off_I(idx,2)-Target.SingleP(idx_struct).TauOffI);
                lowerSD_TInI(iVm) = abs(Tau.Inact_I(idx,2)-Target.SingleP(idx_struct).TauInactI);
                Target.SingleP(idx_struct).TauOnI_SD = (upperSD_TOI(iVm)+lowerSD_TOI(iVm))./2;
                Target.SingleP(idx_struct).TauOffI_SD = (upperSD_TOfI(iVm)+lowerSD_TOfI(iVm))./2;
                Target.SingleP(idx_struct).TauInactI_SD = (upperSD_TInI(iVm)+lowerSD_TInI(iVm))./2;

            end
        end
    end
end

% Tau funcitons wrt V
Vm_val_TV = [-80,-70,-60,-50,-40,-30,-20,-10];
Il_val_TV = fliplr([340,690,2390,5500]);

if any([size(Tau.On_V,1),size(Tau.Inact_V,1)]-(length(Vm_val_TV)*3*length(Il_val_TV)-12)) 
    error('lengts not correct')
end
idx = 0;
for iIl = 1:length(Il_val_TV)
    if Il_val_TV(iIl) == 340
        Vm_val_TV = [-80,-60,-40,-20];
    else
        Vm_val_TV = [-80,-70,-60,-50,-40,-30,-20,-10];
    end

    upperSD_TOV = [];
    lowerSD_TOV = [];
    upperSD_TOfV = [];
    lowerSD_TOfV = [];
    upperSD_TInV = [];
    lowerSD_TInV = [];
    for imsd = 1:3

        for iVm = 1:length(Vm_val_TV)
            idx_struct = [Target.SingleP(:).Vm]==Vm_val_TV(iVm)&[Target.SingleP(:).Il]==Il_val_TV(iIl);
            idx = idx+1;
            if imsd==1
                if ~any(idx_struct)
                    Target.SingleP(end+1).TauOnV = Tau.On_V(idx,2);
                    Target.SingleP(end).TauInactV = Tau.Inact_V(idx,2);
                    Target.SingleP(end).Vm = Vm_val_TV(iVm);
                    Target.SingleP(end).Il = Il_val_TV(iIl);
                    Target.SingleP(end).OSpstart = 0.01;
                    Target.SingleP(end).OSpd = 0.5;
                    Target.SingleP(end).nsamples = 7;
                else
                    Target.SingleP(idx_struct).TauOnV = Tau.On_V(idx,2);
                    Target.SingleP(idx_struct).TauInactV = Tau.Inact_V(idx,2);
                end
            elseif imsd==2

                upperSD_TOV(iVm) = abs(Tau.On_V(idx,2)-Target.SingleP(idx_struct).TauOnV);
                upperSD_TInV(iVm) = abs(Tau.Inact_V(idx,2)-Target.SingleP(idx_struct).TauInactV);

            else
                lowerSD_TOV(iVm) = abs(Tau.On_V(idx,2)-Target.SingleP(idx_struct).TauOnV);
                lowerSD_TInV(iVm) = abs(Tau.Inact_V(idx,2)-Target.SingleP(idx_struct).TauInactV);
                Target.SingleP(idx_struct).TauOnV_SD = (upperSD_TOV(iVm)+lowerSD_TOV(iVm))./2;
                Target.SingleP(idx_struct).TauInactV_SD = (upperSD_TInV(iVm)+lowerSD_TInV(iVm))./2;
            end
        end
    end
end

% Tau funcitons wrt V
Vm_val_TV = [-80,-60,-40,-20];
Il_val_TV = fliplr([340,690,2390,5500]);

if any(size(Tau.Off_V,1)-(length(Vm_val_TV)*3*length(Il_val_TV)))
    error('lengts not correct')
end
idx = 0;
for iIl = 1:length(Il_val_TV)
    upperSD_TOV = [];
    lowerSD_TOV = [];
    upperSD_TOfV = [];
    lowerSD_TOfV = [];
    upperSD_TInV = [];
    lowerSD_TInV = [];
    for imsd = 1:3

        for iVm = 1:length(Vm_val_TV)
            idx_struct = [Target.SingleP(:).Vm]==Vm_val_TV(iVm)&[Target.SingleP(:).Il]==Il_val_TV(iIl);
            idx = idx+1;
            if imsd==1
                if ~any(idx_struct)
                    Target.SingleP(end+1).TauOffV = Tau.Off_V(idx,2);
                    Target.SingleP(end).Vm = Vm_val_Ir(iVm);
                    Target.SingleP(end).Il = Il_val_Ir(iIl);
                    Target.SingleP(end).OSpstart = 0.01;
                    Target.SingleP(end).OSpd = 0.5;
                    Target.SingleP(end).nsamples = 7;
                else
                    Target.SingleP(idx_struct).TauOffV = Tau.Off_V(idx,2);
                end
            elseif imsd==2

                upperSD_TOfV(iVm) = abs(Tau.Off_V(idx,2)-Target.SingleP(idx_struct).TauOffV);

            else
                lowerSD_TOfV(iVm) = abs(Tau.Off_V(idx,2)-Target.SingleP(idx_struct).TauOffV);
                Target.SingleP(idx_struct).TauOffV_SD = (upperSD_TOfV(iVm)+lowerSD_TOfV(iVm))./2;
            end
        end
    end
end
% combine TauI and Vs
TauNames = {'TauOn','TauOff','TauInact'};
for istruct = 1:length(Target.SingleP)
    for iTn = 1:length(TauNames)
        TV = Target.SingleP(istruct).([TauNames{iTn},'V']);
        TI = Target.SingleP(istruct).([TauNames{iTn},'I']);
        TV_SD = Target.SingleP(istruct).([TauNames{iTn},'V_SD']);
        TI_SD = Target.SingleP(istruct).([TauNames{iTn},'I_SD']);
        if (isempty(TV) || isnan(TV)) && ~(isempty(TI) || isnan(TI))
            Target.SingleP(istruct).(TauNames{iTn}) = TI;
        elseif ~(isempty(TV) || isnan(TV)) && (isempty(TI) || isnan(TI))
            Target.SingleP(istruct).(TauNames{iTn}) = TV;
        elseif ~(isempty(TV) || isnan(TV)) && ~(isempty(TI) || isnan(TI))
            Target.SingleP(istruct).(TauNames{iTn}) = (TV+TI)/2;
        end

        if (isempty(TV_SD) || isnan(TV_SD)) && ~(isempty(TI_SD) || isnan(TI_SD))
            Target.SingleP(istruct).([TauNames{iTn},'_SD']) = TI_SD;
        elseif ~(isempty(TV_SD) || isnan(TV_SD)) && (isempty(TI_SD) || isnan(TI_SD))
            Target.SingleP(istruct).([TauNames{iTn},'_SD']) = TV_SD;
        elseif ~(isempty(TV_SD) || isnan(TV_SD)) && ~(isempty(TI_SD) || isnan(TI_SD))
            Target.SingleP(istruct).([TauNames{iTn},'_SD']) = (TV_SD+TI_SD)/2;
        end

    end
end
Target.SingleP = rmfield(Target.SingleP,{'TauOnI','TauOnV','TauOnI_SD','TauOnV_SD',...
    'TauOffI','TauOffV','TauOffI_SD','TauOffV_SD','TauInactI','TauInactV','TauInactI_SD','TauInactV_SD'});

fnSP = fieldnames(Target.SingleP);
for ifnSP = 1:length(fnSP)
    for iF = 1:length(Target.SingleP)
        if isempty(Target.SingleP(iF).(fnSP{ifnSP}))
            Target.SingleP(iF).(fnSP{ifnSP}) = nan;
        end
    end
end
%% TauRecov
Tau.Recov_I=[439.2593434, 2.0000000
    1601.7681291, 1.8959108
    3435.8725538, 2.0297398
    443.4584095, 2.9814126
    1602.0903942, 2.8921933
    3440.1798432, 3.3457249
    479.1552807, 5.3382900
    1600.6858956, 4.5501859
    3440.4275544, 4.1115242
    440.9283880, 7.1598513
    1601.6166164, 7.4275093
    3440.5405877, 4.4609665
    ];
Tau.Recov_V=[-79.8828697, 1.9708665
    -60.0000000, 3.3600426
    -40.0000000, 4.0654971
    -20.1464129, 4.4365127
    -79.8535871, 1.8399361
    -60.0000000, 2.8582244
    -40.0000000, 4.5527699
    -20.0878477, 7.3673792
    -79.8243045, 1.9126421
    -60.0000000, 2.9818608
    -40.0000000, 5.3454971
    -20.1171303, 7.1055823
    ];
%% extract Tau REcov
Vm_val_TI = [-80,-60,-40,-20];
Il_val_TI = [440,1600,3440];
intervals = [0.5,3,5,10,15];
tend = (intervals+0.5)*2+0.01;
nrpulses = 2;

if any([size(Tau.Recov_I,1)]-(length(Vm_val_TI)*length(Il_val_TI)))
    error('lengts not correct')
end
idx = 0;
for iVm = 1:length(Vm_val_TI)

    for iIl = 1:length(Il_val_TI)
        idx = idx+1;
        Target.TauRecov(idx).TauRecovI = Tau.Recov_I(idx,2);
        Target.TauRecov(idx).Vm = Vm_val_TI(iVm);
        Target.TauRecov(idx).Il = Il_val_TI(iIl);
        Target.TauRecov(idx).Pstart = 0.01;
        Target.TauRecov(idx).PD = 0.5;
        Target.TauRecov(idx).nsamples = 4;
        Target.TauRecov(idx).Intervals = intervals;
        Target.TauRecov(idx).nrpulses = nrpulses;
        Target.TauRecov(idx).tend = tend;

    end

end
% Tau funcitons wrt V
Vm_val_TV = [-80,-60,-40,-20];
Il_val_TV = fliplr([440,1600,3440]);

if any([size(Tau.Recov_V,1)]-(length(Vm_val_TV)*length(Il_val_TV)))
    error('lengts not correct')
end
idx = 0;
for iIl = 1:length(Il_val_TV)
    for iVm = 1:length(Vm_val_TV)
        idx_struct = [Target.TauRecov(:).Vm]==Vm_val_TV(iVm)&[Target.TauRecov(:).Il]==Il_val_TV(iIl);
        idx = idx+1;

        if ~any(idx_struct)
            Target.TauRecov(end+1).TauRecovV = Tau.Recov_V(idx,2);
            Target.TauRecov(end).Vm = Vm_val_TV(iVm);
            Target.TauRecov(end).Il = Il_val_TV(iIl);
            Target.TauRecov(end).Pstart = 0.01;
            Target.TauRecov(end).PD = 0.5;
            Target.TauRecov(end).nsamples = 4;
            Target.TauRecov(end).Intervals = intervals;
            Target.TauRecov(end).nrpulses = nrpulses;
            Target.TauRecov(end).tend = tend;
        else
            Target.TauRecov(idx_struct).TauRecovV = Tau.Recov_V(idx,2);
        end
    end
end


% combine TauI and Vs
TauNames = {'TauRecov'};
for istruct = 1:length(Target.TauRecov)
    for iTn = 1:length(TauNames)
        TV = Target.TauRecov(istruct).([TauNames{iTn},'V']);
        TI = Target.TauRecov(istruct).([TauNames{iTn},'I']);
        if (isempty(TV) || isnan(TV)) && ~(isempty(TI) || isnan(TI))
            Target.TauRecov(istruct).(TauNames{iTn}) = TI;
        elseif ~(isempty(TV) || isnan(TV)) && (isempty(TI) || isnan(TI))
            Target.TauRecov(istruct).(TauNames{iTn}) = TV;
        elseif ~(isempty(TV) || isnan(TV)) && ~(isempty(TI) || isnan(TI))
            Target.TauRecov(istruct).(TauNames{iTn}) = (TV+TI)/2;
        end
    end
end
Target.TauRecov = rmfield(Target.TauRecov,{'TauRecovI','TauRecovV',});

fnTR = fieldnames(Target.TauRecov);
for ifnTR = 1:length(fnTR)
    for iF = 1:length(Target.TauRecov)
        if isempty(Target.TauRecov(iF).(fnTR{ifnTR}))
            Target.TauRecov(iF).(fnTR{ifnTR}) = nan;
        end
    end
end


if input('save (1/0)')

    save(['TargetH134R_',datestr(now,'yymmdd'),'.mat'],'Target')
end