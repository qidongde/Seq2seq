clear;
load('UHSAS_2024.mat');

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

save( 'Uhsas_2024_hourly.mat', 'sum_median_hourly', 'time_vec_hourly',"size_dist_median_hourly");
