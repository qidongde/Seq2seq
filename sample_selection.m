%% Hourly median calculate
clear;
load('UHSAS_2014.mat');
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

save( 'Uhsas_2014_hourly.mat', 'sum_median_hourly', 'time_vec_hourly',"size_dist_median_hourly","Dp_");

%% uhsas and cpc data match
clear;

load('Uhsas_2014_hourly.mat');
load('CN_hourly_2013_2024.mat');
cpc_time = [time_all(:,1:3),time_all(:,4)-0.5];

% Data match
[logic,index] = ismember(time_vec_hourly, cpc_time, 'rows');
logic = logic(logic ~= 0);
index = index(index ~= 0);
uhsas_data_match = sum_median_hourly(logic,:);
time_label = time_vec_hourly(logic,:);
cpc_data_match = CN_all(index,1);

% Remove sample including NaN
nan_label = isnan(cpc_data_match);
uhsas_data_match = uhsas_data_match(~nan_label);
time_label = time_label(~nan_label,:);
cpc_data_match = cpc_data_match(~nan_label);

save( 'data_match_2014.mat', 'uhsas_data_match', 'cpc_data_match',"time_label");

%% visualization
clear;

load('data_match_2014.mat');
load('Uhsas_2014_hourly.mat');

daterange = [datenum(2014,2,17,4,0,0) datenum(2014,5,1)];



title_string = {'UHSAS\_heated','Total number concentration UHSAS vs CPC'};
fig = figure;
set(fig,'Color','w','Position',[100 100 1200 800])

% UHSAS_heated
ax1 = axes('Position',[0.08 0.08 0.80 0.25]);
time_uhsas = datenum([time_vec_hourly,zeros(size(time_vec_hourly,1),2)])
idx_uhsas = find(time_uhsas>=daterange(1) & time_uhsas<daterange(2)+1/24);

tmp1 = size_dist_median_hourly(idx_uhsas,:)
PC1 = pcolor(time_uhsas(idx_uhsas),Dp_,tmp1');

set(PC1,'EdgeColor','none')
caxis([0 1200])
ax1.YScale = 'log';
set(ax1,'FontSize',12)
xlim(daterange)
datetick('x','keeplimits')
% ax1.XTickLabel = '';
ax1.YLabel.String = 'D_p (nm)';
ax1.YLabel.FontSize = 15;
ax1.XAxis.MinorTick = 'on';
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
ax2 = axes('Position',[0.08 0.5 0.8 0.25]);
time_stamp = datenum([time_label,zeros(size(time_label,1),2)])
idx_label = find(time_stamp>=daterange(1) & time_stamp<daterange(2)+1/24);
plot(time_stamp(idx_label),uhsas_data_match(idx_label)/cpc_data_match(idx_label))
set(ax2,'FontSize',12)
datetick('x')
xlim(daterange)
ax2.XTickLabel = '';
ax2.YLabel.String = 'UHSAS/CPC';
ax2.YLabel.FontSize = 15;
ax2.XAxis.MinorTick = 'on';
ax2.XAxis.MinorTickValues = [daterange(1):30:daterange(2)];
ax2.XAxis.TickDirection = 'out';
ax2.YAxis.MinorTick = 'on';
ax2.YAxis.TickDirection = 'out';
ax2.YLim = [0 1];
title(title_string{2},'FontSize',15)

linkaxes([ax1 ax2],'x')


% savefig('figname')


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
