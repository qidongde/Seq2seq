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

data_label1 = hours(datetime(uhsas_time_vec_hourly)-t1)+1
uhsas_size_dist_full = NaN(numel(t),99)
uhsas_size_dist_full(data_label1,:)=uhsas_size_dist

data_label2 = hours(datetime(cn_time_label_match)-t1)+1
cn_cpc_data_match_full = NaN(numel(t),1)
cn_uhsas_data_match_full = NaN(numel(t),1)
cn_cpc_data_match_full(data_label2,:)=cn_cpc_data_match
cn_uhsas_data_match_full(data_label2,:)=cn_uhsas_data_match

data_label3 = hours(datetime(bse_time_label)-t1)+0.5
bse_neph_data_match_full = NaN(numel(t),1)
bse_uhsas_data_match_full = NaN(numel(t),1)
bse_neph_data_match_full(data_label3,:)=bse_neph_data_match
bse_uhsas_data_match_full(data_label3,:)=bse_uhsas_data_match

save( 'agg_data_2024.mat', 'uhsas_size_dist_full','uhsas_Dp_bins', ...
    'uhsas_time','cn_uhsas_data_match_full', 'cn_cpc_data_match_full', ...
    'bse_uhsas_data_match_full', 'bse_neph_data_match_full');

%% visualization
clear;

load('agg_data_2024.mat');

% Sys error filter
sys_error_label = zeros(size(uhsas_time));

% hyperparameter
a = 0.25
p = [0.25,0.75];

% CN filter
cn_ratio = cn_uhsas_data_match_full./cn_cpc_data_match_full;
q_cn = quantile(cn_ratio,p);
IQR_cn = q_cn(2)-q_cn(1);
lower_bound_cn = q_cn(1) - a*IQR_cn;
upper_bound_cn = q_cn(2) + a*IQR_cn;
CN_label = or(cn_ratio>upper_bound_cn,cn_ratio<lower_bound_cn);

% Bse filter
bse_ratio = bse_uhsas_data_match_full./bse_neph_data_match_full;
q_bse = quantile(bse_ratio,p);
IQR_bse = q_bse(2)-q_bse(1);
lower_bound_bse = q_bse(1) - a*IQR_bse;
upper_bound_bse = q_bse(2) + a*IQR_bse;
Bse_label = or(bse_ratio>upper_bound_bse,bse_ratio<lower_bound_bse);

label_agg = and(CN_label,Bse_label)
sys_error_label(label_agg)=1


daterange = [datenum(2024,1,1,0,0,0) datenum(2024,6,1,0,0,0)];


title_string = {'UHSAS heatmap','CN', 'Bse'};
fig = figure;
set(fig,'Color','w','Position',[100 100 1800 800])



% UHSAS_heated
ax1 = axes('Position',[0.08 0.08 0.80 0.25]);
% ax1 = nexttile;
time_uhsas_heated = datenum(uhsas_time)
idx_uhsas_heated = find(time_uhsas_heated>=daterange(1) & time_uhsas_heated<daterange(2));

sys_label=sys_error_label(idx_uhsas_heated);
x_all = time_uhsas_heated(idx_uhsas_heated)
tmp1 = uhsas_size_dist_full(idx_uhsas_heated,:)
PC1 = pcolor(x_all,uhsas_Dp_bins,tmp1');

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

cn_y = cn_ratio(idx_uhsas_heated)
normal = find(sys_label==0);
sys_error = find(sys_label==1);

plot(x_all(normal),cn_y(normal));
hold on
scatter(x_all(sys_error),cn_y(sys_error),'x');
legend('Normal','Sys error')

set(ax2,'FontSize',12)
datetick('x')
xlim(daterange)
ax2.XTickLabel = '';
ax2.YLabel.String = 'UHSAS/CPC';
ax2.YLabel.FontSize = 15;
ax2.XAxis.MinorTick = 'off';
ax2.XAxis.MinorTickValues = [daterange(1):30:daterange(2)];
ax2.XAxis.TickDirection = 'out';
ax2.YAxis.MinorTick = 'off';
ax2.YAxis.TickDirection = 'out';

title(title_string{2},'FontSize',15)


% 'UHSAS vs Reph'
ax3 = axes('Position',[0.08 0.7 0.80 0.25]);

bse_y = bse_ratio(idx_uhsas_heated)

plot(x_all(normal),bse_y(normal));
hold on
scatter(x_all(sys_error),bse_y(sys_error),'x');
legend('Normal','Sys error')

set(ax2,'FontSize',12)
datetick('x')
xlim(daterange)
ax3.XTickLabel = '';
ax3.YLabel.String = 'UHSAS/Neph';
ax3.YLabel.FontSize = 15;
ax3.XAxis.MinorTick = 'off';
ax3.XAxis.MinorTickValues = [daterange(1):30:daterange(2)];
ax3.XAxis.TickDirection = 'out';
ax3.YAxis.MinorTick = 'off';
ax3.YAxis.TickDirection = 'out';
% ax3.YLim = [0 1.5];
title(title_string{3},'FontSize',15)

linkaxes([ax1 ax2 ax3],'x')

% savefig('figname')

%% test
clear;

load('CN_hourly_2013_2024.mat');

daterange = [datenum(2015,1,1,0,0,0) datenum(2015,8,1,0,0,0)];

time_vec = [time_all,zeros(size(time_all,1),2)];
time_stamp = datenum(time_vec);

idx_label = find(time_stamp>=daterange(1) & time_stamp<daterange(2));

plot(time_stamp(idx_label),CN_all(idx_label))
datetick('x','yyyy-mm-dd')
%% heated update

clear;

load('agg_data_2014.mat');

t1 = datetime(2014,1,1,0,0,0);
t2 = datetime(2015,1,1,0,0,0);
t = t1:hours(1):t2
t=t'
time_vec = datevec(t)
data_label = hours(datetime(uhsas_time_vec_hourly)-t1)+1

uhsas_size_dist_full = NaN(numel(t),99)
uhsas_size_dist_full(data_label,:)=uhsas_size_dist

daterange = [datenum(2014,6,1,0,0,0) datenum(2014,12,1,0,0,0)];

fig = figure;
set(fig,'Color','w','Position',[100 100 1200 800])
title_string = {'UHSAS heated','CN', 'Bse'};

tiledlayout(2,1)


% UHSAS_heated
ax1 = axes('Position',[0.08 0.08 0.80 0.25]);
% ax1 = nexttile;
time_uhsas_heated = datenum(time_vec)
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

%% Filter
clear;

load('agg_data_2014.mat');

sys_error_label = zeros(size(uhsas_time));

% hyperparameter
a = 3
p = [0.25,0.75];

% CN filter
cn_ratio = cn_uhsas_data_match_full./cn_cpc_data_match_full;
q_cn = quantile(cn_ratio,p);
IQR_cn = q_cn(2)-q_cn(1);
lower_bound_cn = q_cn(1) - a*IQR_cn;
upper_bound_cn = q_cn(2) + a*IQR_cn;
CN_label = or(cn_ratio>upper_bound_cn,cn_ratio<lower_bound_cn);

% Bse filter
bse_ratio = bse_uhsas_data_match_full./bse_neph_data_match_full;
q_bse = quantile(bse_ratio,p);
IQR_bse = q_bse(2)-q_bse(1);
lower_bound_bse = q_bse(1) - a*IQR_bse;
upper_bound_bse = q_bse(2) + a*IQR_bse;
Bse_label = or(bse_ratio>upper_bound_bse,bse_ratio<lower_bound_bse);

label_agg = and(CN_label,Bse_label)
sys_error_label(label_agg)=1

linkaxes([ax1 ax2 ax3],'x')

% savefig('figname')
