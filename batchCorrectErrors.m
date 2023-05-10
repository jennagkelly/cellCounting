function batchCorrectErrors

channels={'GABA','Kv3','PV'};

wkspctofix='m652b542n99_V2_DKPG_ser1_12162019.mat';
applyManCorrectionsa(wkspctofix,channels)
wkspctofix='m652b542n99_V2_DKPG_ser1_12162019CORR.mat';
applyManCorrectionsb(wkspctofix,channels)


wkspctofix='m652b552n103_V2_DKPG_ser1_11132019.mat';
applyManCorrectionsa(wkspctofix,channels)
wkspctofix='m652b552n103_V2_DKPG_ser1_11132019CORR.mat';
applyManCorrectionsb(wkspctofix,channels)


wkspctofix='m652b552n103_V2_DKPG_ser2_11132019.mat';
applyManCorrectionsa(wkspctofix,channels)
wkspctofix='m652b552n103_V2_DKPG_ser2_11132019CORR.mat';
applyManCorrectionsb(wkspctofix,channels)


wkspctofix='M652b542n148_V2_DKPG_ser1_11102019.mat';
applyManCorrectionsa(wkspctofix,channels)
wkspctofix='M652b542n148_V2_DKPG_ser1_11102019CORR.mat';
applyManCorrectionsb(wkspctofix,channels)


wkspctofix='m652b542n148_V2_DKPG_ser2_12192019.mat';
applyManCorrectionsa(wkspctofix,channels)
wkspctofix='m652b542n148_V2_DKPG_ser2_12192019CORR.mat';
applyManCorrectionsb(wkspctofix,channels)


wkspctofix='M652b552n160_V2_DKPG_ser1_12162019.mat';
applyManCorrectionsa(wkspctofix,channels)
wkspctofix='M652b552n160_V2_DKPG_ser1_12162019CORR.mat';
applyManCorrectionsb(wkspctofix,channels)


wkspctofix='M652b552n160_V2_DKPG_ser2_12162019.mat';
applyManCorrectionsa(wkspctofix,channels)
wkspctofix='M652b552n160_V2_DKPG_ser2_12162019CORR.mat';
applyManCorrectionsb(wkspctofix,channels)


