%% Hourly median calculate
clear;
load('UHSAS_2023.mat');
Dp_ = mean(Dp_bounds,1)'

% findgroups
time_vec = datevec(Time_UHSAS);
[G,year,month,day,hour] = findgroups(time_vec(:,1),time_vec(:,2),time_vec(:,3),time_vec(:,4));
time_vec_hourly = [year month day hour];

% Hourly total concentration median
Dp_range = log(Dp_bounds(2,:))-log(Dp_bounds(1,:));
Dp_sum = dN_dlogDp * Dp_range';
sum_median_hourly = splitapply(@median,Dp_sum,G);

% Hourly size distribution median
size_dist_median_hourly = []
for i=1:99
    tmp = splitapply(@median,dN_dlogDp(:,i),G);
    size_dist_median_hourly = [size_dist_median_hourly,tmp]
end

save( 'Uhsas_2023_hourly.mat', 'sum_median_hourly', 'time_vec_hourly',"size_dist_median_hourly","Dp_");

%% uhsas and cpc data match
clear;

load('Uhsas_2024_hourly.mat');
load('CN_hourly_2013_2024.mat');
cpc_time = [time_all(:,1:3),time_all(:,4)-0.5];

% Data match
[logic,index] = ismember(time_vec_hourly, cpc_time, 'rows');
% logic = find(logic);
index = index(index ~= 0);
uhsas_data_match = sum_median_hourly(logic,:);
time_label = time_vec_hourly(logic,:);
cpc_data_match = CN_all(index,1);

% Remove sample including NaN
nan_label = isnan(cpc_data_match);
cn_uhsas_data_match = uhsas_data_match(~nan_label);
cn_cpc_data_match = cpc_data_match(~nan_label);
time_label = time_label(~nan_label,:);
cn_time_label_match = [time_label,zeros(size(time_label,1),2)]



save( 'cn_data_match_2024.mat', 'cn_uhsas_data_match', 'cn_cpc_data_match',"cn_time_label_match");
%% uhsas vs cpc test

time_stamp = datenum(cn_time_label_match);
time_stamp2 = datenum([cpc_time,zeros(size(cpc_time,1),2)]);
daterange = [datenum(2018,4,1,4,0,0) datenum(2018,9,30)];
idx_label = find(time_stamp>=daterange(1) & time_stamp<daterange(2)+1/24);
idx_label2 = find(time_stamp2>=daterange(1) & time_stamp2<daterange(2)+1/24);

tiledlayout(2,1)

ax1 = nexttile;
scatter(ax1,time_stamp2(idx_label2),CN_all(idx_label2,1))
xlim(daterange)
datetick('x','keeplimits')
ax1.XAxis.MinorTick = 'off';
ax1.XAxis.MinorTickValues = [daterange(1):30:daterange(2)];
ax1.XAxis.TickDirection = 'out';
title('Raw','FontSize',15)

ax2 = nexttile;
scatter(ax2,time_stamp(idx_label),cn_cpc_data_match(idx_label),'filled','d')
ax2.XAxis.MinorTick = 'off';
ax2.XAxis.MinorTickValues = [daterange(1):30:daterange(2)];
ax2.XAxis.TickDirection = 'out';
xlim(daterange)
datetick('x','keeplimits')
title('Match','FontSize',15)


%% uhsas and Reph data match
clear;

load('Reph_data_2024.mat');
load('uhsas_data_2024.mat');


% Data match
[logic,index] = ismember(time_vec_hourly_, Neph_time_vec, 'rows');
% logic = find(logic);
index = index(index ~= 0);
uhsas_data_match = Bsp_Mie(logic,:);
time_label = time_vec_hourly_(logic,:);
reph_data_match = Bsp_median_hourly(index);

% Remove sample including NaN
nan_label = isnan(reph_data_match);
bse_uhsas_data_match = uhsas_data_match(~nan_label);
bse_neph_data_match = reph_data_match(~nan_label);
bse_time_label = time_label(~nan_label,:);

save( 'bse_data_match_2024.mat', 'bse_uhsas_data_match', 'bse_neph_data_match',"bse_time_label");

%% Data organize
clear;

load('cn_data_match_2024.mat');
load('bse_data_match_2024.mat');
load('Uhsas_2024_hourly.mat');

uhsas_size_dist = size_dist_median_hourly;
uhsas_Dp_bins = Dp_;
uhsas_time_vec_hourly = [time_vec_hourly,zeros(size(time_vec_hourly,1),2)];

t1 = datetime(2024,1,1,0,0,0);
t2 = datetime(2025,1,1,0,0,0);
t = t1:hours(1):t2
uhsas_time=t'
data_label = hours(datetime(uhsas_time_vec_hourly)-t1)+1

uhsas_size_dist_full = NaN(numel(t),99)
uhsas_size_dist_full(data_label,:)=uhsas_size_dist


save( 'agg_data_2024.mat', 'uhsas_size_dist_full','uhsas_Dp_bins','uhsas_time', ...
    'cn_uhsas_data_match', 'cn_cpc_data_match',"cn_time_label_match", ...
    'bse_uhsas_data_match', 'bse_neph_data_match',"bse_time_label");
    
%% scattering coefficient

clear;

load('Uhsas_2024_hourly.mat');
load('UHSAS_2014.mat','Dp_bounds');

time_vec_hourly_ = [time_vec_hourly,30*ones(size(time_vec_hourly,1),1),zeros(size(time_vec_hourly,1),1)];
time_stamp = datenum(time_vec_hourly_);

% Dp, geometric mean value of boundaries, unit m, (1,99)
Dp = geomean(Dp_bounds,1)/1e9;
% dlnDp, logarithmic difference of particle diameter boundaries, (1,99)
dlnDp = log(Dp_bounds(2,:)./Dp_bounds(1,:));
% n, size distribution (dN/dlnDp) cm^-3, 2 dimensional, (sample_num,99)
n = size_dist_median_hourly;
% m, refractive index, related to components, default=1.5
m = 1.5;
% wl, wavelength, unit m, default=550*1e-9, green light
wl = 550*1e-9;

% scattering coefficient by Mie theory
Bsp_Mie = aerosol_optical_properties(Dp,dlnDp,n,m,wl)*1e6
Bsp_Mie = Bsp_Mie'

save( 'uhsas_data_2024.mat', 'size_dist_median_hourly','sum_median_hourly', "Bsp_Mie",'time_vec_hourly_','time_stamp','Dp_bounds');

%% Neph hourly data
clear;
filename = 'E:\wustl\size distribution\target clean\scattering coe\Neph_2014.mat';
load(filename);

% findgroups
time_vec = datevec(Time_Neph);
[G,year,month,day,hour] = findgroups(time_vec(:,1),time_vec(:,2),time_vec(:,3),time_vec(:,4));
time_vec_hourly = [year month day hour];
Neph_time_vec = [time_vec_hourly,30*ones(size(time_vec_hourly,1),1),zeros(size(time_vec_hourly,1),1)];
Neph_time_stamp = datenum(Neph_time_vec);

% Hourly scattering coef
Bsp_median_hourly = splitapply(@median,Bs_RGB(:,2),G);

save( 'Reph_data_2014.mat', "Bsp_median_hourly",'Neph_time_vec','Neph_time_stamp');

%% compare Bsp from Neph and Mie theory

clear;

load('uhsas_data_2014.mat');
load('Reph_data_2014.mat');

daterange = [datenum(2014,3,17,4,0,0) datenum(2014,4,17,4,0,0)];

ax2 = axes();
% ax2 = axes('Position',[0.08 0.5 0.8 0.25]);

Neph_label = find(Neph_time_stamp>=daterange(1) & Neph_time_stamp<daterange(2)+1/24);
plot(Neph_time_stamp(Neph_label),Bsp_median_hourly(Neph_label));

hold on;
Mie_label = find(time_stamp>=daterange(1) & time_stamp<daterange(2)+1/24);
scatter(time_stamp(Mie_label),Bsp_Mie(Mie_label),'.');
hold off;

legend('Neph','Mie')

set(ax2,'FontSize',12)
datetick('x','yyyy-mm-dd')
% xlim(daterange)
% ax2.XTickLabel = '';
% ax2.YLabel.String = 'Bsp(1/Mm)';
% ax2.YLabel.FontSize = 15;
ax2.XAxis.MinorTick = 'on';
ax2.XAxis.MinorTickValues = [daterange(1):5:daterange(2)];
% ax2.XAxis.TickDirection = 'out';
% ax2.YAxis.MinorTick = 'on';
% ax2.YAxis.TickDirection = 'out';
% ax2.YLim = [0 1];
title('Scattering Coef','FontSize',15)

% linkaxes([ax1 ax2],'x')


%% visualization
clear;

load('agg_data_2014.mat');

daterange = [datenum(2014,1,1,0,0,0) datenum(2014,6,1,0,0,0)];


title_string = {'UHSAS heated','CN', 'Bse'};
fig = figure;
set(fig,'Color','w','Position',[100 100 1200 800])

tiledlayout(2,1)


% UHSAS_heated
ax1 = axes('Position',[0.08 0.08 0.80 0.25]);
% ax1 = nexttile;
time_uhsas_heated = datenum(uhsas_time)
idx_uhsas_heated = find(time_uhsas_heated>=daterange(1) & time_uhsas_heated<daterange(2));

tmp1 = uhsas_size_dist_full(idx_uhsas_heated,:)
PC1 = pcolor(time_uhsas_heated(idx_uhsas_heated),uhsas_Dp_bins,tmp1');

set(PC1,'EdgeColor','none')
caxis([0 1200])
ax1.YScale = 'log';
set(ax1,'FontSize',12)
xlim(daterange)
datetick('x','yyyy-mm-dd')
% ax1.XTickLabel = '';
ax1.YLabel.String = 'D_p (nm)';
ax1.YLabel.FontSize = 15;
ax1.XAxis.MinorTick = 'off';
ax1.XAxis.MinorTickValues = [daterange(1):30:daterange(2)];
ax1.XAxis.TickDirection = 'out';
ax1.YAxis.MinorTick = 'on';
ax1.YAxis.TickDirection = 'out';
ax1.YLim = [60 1000];
title(title_string{1},'FontSize',15)

h1 = colorbar('v');
h1.Limits = [0 1200];
colormap(jet)
h1.Position = [0.89 0.08 0.015 0.25];
h1.FontSize = 12;
h1.Label.String = 'd{\itN}/d{\itlogD_p} (cm^{–3})';
h1.Label.FontSize = 15;
h1.Ticks = [0:300:1200];
h1.TickDirection = 'out';
h1.TickLength = 0.03;

% 'UHSAS vs CPC'
ax2 = axes('Position',[0.08 0.38 0.80 0.25]);
% ax2 = nexttile;
time_cn = datenum(cn_time_label_match);
idx_label_cn = find(time_cn>=daterange(1) & time_cn<daterange(2));
plot(time_cn(idx_label_cn),cn_uhsas_data_match(idx_label_cn)/cn_cpc_data_match(idx_label_cn))
set(ax2,'FontSize',12)
datetick('x')
xlim(daterange)
ax2.XTickLabel = '';
ax2.YLabel.String = 'UHSAS/CPC';
ax2.YLabel.FontSize = 15;
ax2.XAxis.MinorTick = 'off';
ax2.XAxis.MinorTickValues = [daterange(1):30:daterange(2)];
ax2.XAxis.TickDirection = 'out';
ax2.YAxis.MinorTick = 'on';
ax2.YAxis.TickDirection = 'out';
% ax2.YLim = [0 0.8];
title(title_string{2},'FontSize',15)


% 'UHSAS vs Reph'
ax3 = axes('Position',[0.08 0.7 0.80 0.25]);
% ax3 = nexttile;
time_stamp_bse = datenum(bse_time_label)
idx_label_bse = find(time_stamp_bse>=daterange(1) & time_stamp_bse<daterange(2));
plot(time_stamp_bse(idx_label_bse),bse_uhsas_data_match(idx_label_bse)/bse_neph_data_match(idx_label_bse))
set(ax2,'FontSize',12)
datetick('x')
xlim(daterange)
ax3.XTickLabel = '';
ax3.YLabel.String = 'UHSAS/Neph';
ax3.YLabel.FontSize = 15;
ax3.XAxis.MinorTick = 'off';
ax3.XAxis.MinorTickValues = [daterange(1):30:daterange(2)];
ax3.XAxis.TickDirection = 'out';
ax3.YAxis.MinorTick = 'on';
ax3.YAxis.TickDirection = 'out';
% ax3.YLim = [0 1.5];
title(title_string{3},'FontSize',15)

linkaxes([ax1 ax2 ax3],'x')

% savefig('figname')
