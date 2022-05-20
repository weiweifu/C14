dir0 = '../Mask_TTM/';
fname01 = 'MSK_240km_dst.mat';
fname02 = 'MSK_60km_src.mat';

load ../Mask_TTM/Lump_Spray.mat S;
load ../Mask_TTM/Lump_Spray2.mat L;

grd1 = load([dir0 fname01]);
grd2 = load([dir0 fname02]);

aa1 = grd1.GRD.TEMP';

aa2 = grd2.GRD.TEMP';

% figure 1 low-res sst
figure(1);
plotaa(grd1.GRD,aa1,1);
title('SST (QU240km)');
export_fig('SST_c','-png');

figure(2);
plotaa(grd2.GRD,aa2,1);
title('SST (oEC60to30)');
export_fig('SST_f','-png');

bb1 = grd1.GRD.TEMP';
bb1(grd1.GRD.iocn) = L*aa2(grd2.GRD.iocn);

figure(3)
plotaa(grd1.GRD,bb1,1);
title('SST (oEC60to30 --> QU240km)');
export_fig('SST_f-c','-png');

% figure 3 low-res remap sst

bb2 = grd2.GRD.TEMP';
bb2(:) = 0;
bb2(grd2.GRD.iocn) = S*aa1(grd1.GRD.iocn);
figure(4)
plotaa(grd2.GRD,bb2,1);
title('SST (QU240km --> oEC60to30)');
export_fig('SST_c-f','-png');
