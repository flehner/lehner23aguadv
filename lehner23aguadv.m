% ---------------------------
% lehner23aguadv.m
% ---------------------------
% Code to produce Fig. 1 in Lehner et al. (2023)
%
% Citation:
% Lehner, F., E. Hawkins, R. Sutton, A. G. Pendergrass, F. C. Moore (2023):
% New potential to reduce uncertainty in regional climate projections by combining physical and socio-economic constraints
% AGU Advances, DOI: XXX
%
% Notes on code:
% - Only does calculations and plotting, not pre-processing
% - Requires pre-processed monthly mean data, specifically:
%   - Connected CMIP6 simulations from historical, ssp126, ssp245, ssp370, ssp585, bilinearly regridded to 2.5°x2.5°
%   - Raw CMIP6 data: https://esgf-node.llnl.gov/search/cmip6/
% ---------------------------

close all
clear all

% -- Main paramaters to edit ---------------------------------------------------

% -- Path to data:
pathin    = '~/Dropbox/work/';

% -- 1st variable (predictor) options:
vari      = 'tas';
% seasons     = {'DJF','JJA','annual','JJAS'};
seas      = 'annual';

% -- 2nd variable (predictand) options:
vari2     = 'tas'; % pr tas siconc
seas2     = 'JJA'; % annual JJA S

% -- region (search for "regions" for full list):
r         = 30;

% -- land-only data or global?
lnd_only  = 0; % 1=yes, 0=no ; for most cases it's pre-defined

% ------------------------------------------------------------------------------

% -- emissions scenarios:
scen_cmip6   = {'ssp119','ssp126','ssp245','ssp370','ssp585'};
% -- CMIP6 models:
models_cmip6 = {'ACCESS-CM2','ACCESS-ESM1-5','AWI-CM-1-1-MR','BCC-CSM2-MR','CAMS-CSM1-0',...
                'CESM2','CESM2-WACCM','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1',...
                'CanESM5','CanESM5-CanOE','EC-Earth3','EC-Earth3-Veg','FGOALS-f3-L',...
                'FGOALS-g3','FIO-ESM-2-0','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G',...
                'HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G',...
                'MCM-UA-1-0','MIROC-ES2L','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR',...
                'MRI-ESM2-0','NESM3','NorESM2-LM','NorESM2-MM','UKESM1-0-LL', 'EC-Earth3-Veg-LR'};
% -- excluded CMIP6 models:
models_cmip6(ismember(models_cmip6,'NorESM2-LM')) = []; % pr: ssp585 r1i1p1f1 exists, but the historical r1i1p1f1 only starts in 1950 (?)
models_cmip6(ismember(models_cmip6,'FIO-ESM-2-0')) = []; % no ssp370
models_cmip6(ismember(models_cmip6,'HadGEM3-GC31-LL')) = []; % no ssp370
models_cmip6(ismember(models_cmip6,'NESM3')) = []; % no ssp370
models_cmip6(ismember(models_cmip6,'KACE-1-0-G')) = []; % piControl only 150 years

if length(scen_cmip6)==5
  % -- the following models have ssp119:
  models_cmip6 = {'CAMS-CSM1-0','CNRM-ESM2-1','CanESM5','EC-Earth3-Veg-LR','EC-Earth3-Veg','FGOALS-g3','GFDL-ESM4','GISS-E2-1-G','GISS-E2-1-H','IPSL-CM6A-LR','MIROC-ES2L','MIROC6','MPI-ESM1-2-LR','MRI-ESM2-0','UKESM1-0-LL'};
  if strcmp(vari2,'siconc')==1
    % -- and if you want ssp119 and siconc as 2nd variable, then it's these:
    models_cmip6 = {'CAMS-CSM1-0','CNRM-ESM2-1','CanESM5','EC-Earth3-Veg','FGOALS-g3','GFDL-ESM4','GISS-E2-1-H','IPSL-CM6A-LR','MIROC-ES2L','MIROC6','MPI-ESM1-2-LR','MRI-ESM2-0','UKESM1-0-LL'};
  end
end

ensmem_cmip6  = ones(length(models_cmip6),1);

% -- start and end year of CMIP6 simulations:
start_cmip6   = 1850;
ende_cmip6    = 2100;
% -- common time period:
start         = 1950;
ende          = 2099;
time          = start:ende;
% -- reference period:
refstart      = 2001; % 2001 1995 1986 1971 1950 2001 1979 1991
refende       = 2020; % 2020 2014 2005 2000 1979 2030 2000 2020
reflength     = refende-refstart+1;
% -- runnig mean length in years:
wl            = 10; % 1 10
% -- grid:
lati          = [-88.75:2.5:88.75];
loni          = [1.25:2.5:358.75];

% -- area and land mask:
filein    = [pathin 'cmip5-ng/area_g025.nc'];
area      = ncread([filein],'AREA');
filein    = [pathin 'cmip5-ng/landfrac_g025.nc'];
land      = ncread([filein],'LANDFRAC');
id        = isnan(land);
land(id)  = 0;
landarea  = land.*area;

% -- regions:
regions     = {'global','uk','sahel','europe','southern_ocean',...
               'se_asia','us_southwest','nino34','north_america','kolkata',...
               'dallas','alaska','seattle','upper_colorado','sydney',...
               'NH','SH','arctic','AA','zurich',...
               'southern_europe','northern_europe','denver','anchorage','india',...
               'sahara','NoCal','SoCal','Mediterranean','Central_Europe',...
               'Northern_Europe','NCA_nw','UCRB','ireland','NCA_CONUS',...
               'Amazon','Hudson_Bay'};

% -- SREX region definition:
tmp = ncread([pathin 'cmip5-ng/ipcc_srex_regions.nc'],'regions_bin');
% -- longitude flip:
srex_regions = tmp;
srex_regions(1:72,:,:) = tmp(73:end,:,:);
srex_regions(73:end,:,:) = tmp(1:72,:,:);


region = regions{r};
clear('lon2','var0','seas0')
% -- region definitions:
if strcmp(region,'NCA')==1
  lat         = [20 70];
  lon         = [-125 -65]+360;
end
if strcmp(region,'NCA_CONUS')==1
  lat         = [26.1 49.1];
  lon         = [-124.9+360 -67.1+360];
end
if strcmp(region,'Hudson_Bay')==1
  lat         = [55.4 62.8];
  lon         = [-94.5+360 -77.2+360];
end
if strcmp(region,'Amazon')==1
  lat         = [-7.1 4.3];
  lon         = [-69.9+360 -48.3+360];
end
if strcmp(region,'upper_colorado')==1
  lat         = [38.75 41.25];
  lon         = [248.75 251.25 253.75];
    seas0     = 3;
    var0      = 2;
end
if strcmp(region,'UCRB')==1
  lat         = [37 42];
  lon         = [-111 -105]+360;
end
if strcmp(region,'us_southwest')==1
  lat         = [32.5 42.1];
  lon         = [-125 -103]+360;
end
if strcmp(region,'NoCal')==1
  lat         = [37.1 42.1];
  lon         = [-125 -113]+360;
end
if strcmp(region,'SoCal')==1
  lat         = [32.5 37.1];
  lon         = [-125 -113]+360;
end
if strcmp(region,'europe')==1
  lat         = [36.4 70.1];
  lon         = [0 23.4];
  lon2        = [350 360];
end
if strcmp(region,'southern_europe')==1
  lat         = [36.4 50.1];
  lon         = [0 23.4];
  lon2        = [350 360];
end
if strcmp(region,'northern_europe')==1
  lat         = [50.1 70.1];
  lon         = [0 23.4];
  lon2        = [350 360];
end
if strcmp(region,'se_asia')==1
  lat         = [7.9 30.1];
  lon         = [93.2 122.1];
end
if strcmp(region,'india')==1
  lat         = [6.2 30.1];
  lon         = [69.9 87.3];
    seas0     = 2;
    var0      = 2;
end
if strcmp(region,'nino34')==1
  lat         = [-5 5];
  lon         = [-170 -120]+360;
    seas0        = 1;
    var0        = 1;
end
if strcmp(region,'uk')==1
  lat         = [50.2:2.5:58.8];
  lon         = [-10.9:2.5:0]+360;
    seas0       = 3;
    var0        = 1;
end
if strcmp(region,'ireland')==1
  lat         = [52.2:2.5:55.8];
  lon         = [-10.9:2.5:-5]+360;
end
if strcmp(region,'sydney')==1
  lat         = -33.8;
  lon         = 150.0;
end
if strcmp(region,'dallas')==1
  lat         = 32.7;
  lon         = -96.8+360;
end
if strcmp(region,'seattle')==1
  lat         = 49.2;
  lon         = -121.4+360;
    var0      = 2;
    seas0      = 1;
end
if strcmp(region,'sahara')==1
  lat         = 22.5;
  lon         = 9.6;
    var0      = 2;
    seas0     = 3;
end
if strcmp(region,'denver')==1
  lat         = 39.9;
  lon         = -100.1+360;
end
if strcmp(region,'anchorage')==1
  lat         = 61.2;
  lon         = -149.1+360;
end
if strcmp(region,'north_america')==1
  lat         = [22.1 70.5];
  lon         = [-124.4+360 -53.7+360];
end
if strcmp(region,'kolkata')==1
  lat         = 22.6;
  lon         = 88.4;
end
if strcmp(region,'nao')==1
  lat         = 50; % dummy, actual NAO index is loaded later
  lon         = 50; % dummy, actual NAO index is loaded later
end
if strcmp(region,'southern_ocean')==1
  lat         = [-70 -60];
  lon         = [0 360];
    seas0      = 3;
    var0      = 1;
end
if strcmp(region,'sahel')==1
  lat         = [10 20];
  lon         = [0 40];
    var0      = 2;
    seas0     = 2;
end
if strcmp(region,'NH')==1
  lat         = [0 90];
  lon         = [0 360];
end
if strcmp(region,'SH')==1
  lat         = [-90 0];
  lon         = [0 360];
end
if strcmp(region,'alaska')==1
  lat         = [51.6 73.9];
  lon         = [191.2 255.9];
end
if strcmp(region,'zurich')==1
  lat         = 47.5;
  lon         = 8.5;
end
if strcmp(region,'arctic')==1
  lat         = [70 90];
  lon         = [0 360];
end
if strcmp(region,'AA')==1
  lat1        = [0 90];
  lat2        = [75 90];
  lon1        = [0 360];
  lon2        = [0 360];
end
% --
% -- find grid cells to plot in models
clear('tmp','tmp0')
if strcmp(region,'global')==1
  ii      = 1:length(lati);
  jj      = 1:length(loni);
  weights = area(jj,ii);
  weights = weights/nansum(nansum(weights));
  weights_lnd = landarea(jj,ii);
  weights_lnd = weights_lnd/nansum(nansum(weights_lnd));
else if strcmp(region,'Mediterranean')==1
  ii      = 1:length(lati);
  jj      = 1:length(loni);
  weights = area .* squeeze(srex_regions(:,:,12));
  weights = weights/nansum(nansum(weights));
else if strcmp(region,'Central_Europe')==1
  ii      = 1:length(lati);
  jj      = 1:length(loni);
  weights = area .* squeeze(srex_regions(:,:,6));
  weights = weights/nansum(nansum(weights));
else if strcmp(region,'Northern_Europe')==1
  ii      = 1:length(lati);
  jj      = 1:length(loni);
  weights = area .* squeeze(srex_regions(:,:,16));
  weights = weights/nansum(nansum(weights));
else if strcmp(region,'NCA_nw')==1
  % -- NCA regions
  filein  = [pathin 'NCA_masks_g025.nc'];
  nca_region = ncread([filein],'NW');
  weights = area .* nca_region;
  weights = weights/nansum(nansum(weights));
else if strcmp(region,'AA')==1
  for x = 1:length(lat1)
    tmp0    = find(abs(lati-lat1(x))==min(abs(lati-lat1(x))));
    tmp1(x)  = tmp0(1);
  end
  for x = 1:length(lat2)
    tmp0    = find(abs(lati-lat2(x))==min(abs(lati-lat2(x))));
    tmp2(x)  = tmp0(1);
  end
  ii1 = tmp1(1):tmp1(end);
  ii2 = tmp2(1):tmp2(end);
  for x = 1:length(lon1)
    tmp0    = find(abs(loni-lon1(x))==min(abs(loni-lon1(x))));
    tmp(x)  = tmp0(1);
  end
  jj1 = tmp(1):tmp(end);
  jj2 = jj1;
  % -- cut out area --
  weights1 = area(jj1,ii1);
  weights2 = area(jj2,ii2);
  weights1 = weights1/nansum(nansum(weights1));
  weights2 = weights2/nansum(nansum(weights2));
  clear('tmp1','tmp2')
else
  for x = 1:length(lat)
    tmp0    = find(abs(lati-lat(x))==min(abs(lati-lat(x))));
    tmp(x)  = tmp0(1);
  end
  ii = tmp(1):tmp(end);
  for x = 1:length(lon)
    tmp0    = find(abs(loni-lon(x))==min(abs(loni-lon(x))));
    tmp(x)  = tmp0(1);
  end
  jj = tmp(1):tmp(end);
  if exist('lon2')==1
    clear('tmp')
    for x = 1:length(lon2)
      tmp0    = find(abs(loni-lon2(x))==min(abs(loni-lon2(x))));
      tmp(x)  = tmp0(1);
    end
    jj = [jj tmp(1):tmp(end)];
  end
  % -- cut out area --
  if strcmp(region,'nino34')==1 ||...
     strcmp(region,'southern_ocean')==1 ||...
     strcmp(region,'arctic')==1
    weights = area(jj,ii);
    weights = weights/nansum(nansum(weights));
  else
    weights = landarea(jj,ii);
    weights = weights/nansum(nansum(weights));
  end
end
end
end
end
end
end
clear('tmp','tmp0')

% -- print to screen the current parameters
'------------------------------------------------'
['Region: ' region ]
['1st variable: ' vari ]
['1st variable season: ' seas ]
['2nd variable: ' vari2 ]
['2nd variable season: ' seas2 ]
'------------------------------------------------'


% -- Load data --
% -- 1st variable:
if strcmp(vari,'tas')==1
  units       = 'K';
  f           = 1;
  % units       = '\circF';
  % f           = 1.8;
  var_name    = 'Temperature'
  var_name2   = 'temperature';
elseif strcmp(vari,'pr')==1
  units       = '%';
  f           = 86400;
  var_name    = 'Precipitation'
  var_name2   = 'precipitation';
elseif strcmp(vari,'psl')==1
  units       = 'unitless';
  f           = 1;
  var_name    = 'Sea level pressure'
end
if refende <= 2020
  % -- load global mean observations
  obs_raw = squeeze(ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1.185001-202012_ts.nc'],'temperature'));
  obs_raw_am = am(obs_raw);
  if strcmp(vari,'tas')==1 && strcmp(region,'global')==1
    obs_name = 'BEST';
    time_obs = 1850:2022;
    % obs_raw = squeeze(ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1.185001-202012_ts.nc'],'temperature'));
    obs_raw = squeeze(ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1_g025_ts.nc'],'temperature'));
    obs_raw_am = am(obs_raw);
    obs = rm(am(obs_raw),wl);
    obs = obs-nanmean(obs(refstart-1850+1:refende-1850+1));
  else
    obs_name = 'GPCP';
    time_obs = 1986:2020;
    obs_raw  = squeeze(ncread([pathin 'observations/precipitation/gpcp/precip.mon.mean.197901-202012_ts.nc'],'precip'));
    obs_raw  = obs_raw(85:end); % limit to 1986 onwards, based on Angie's assessment
    obs  = rm(am(obs_raw),wl);
    ref  = nanmean(obs(refstart-1986+1:refende-1986+1));
    obs  = ((obs-ref)/ref)*100;
  end
  tmp1                = rm(am(obs_raw),wl)'; % do smoothing
  if wl > 1
    tmp1 = tmp1';
  end
  idx = ~isnan(tmp1);
  fit                 = polyval(polyfit(time_obs(idx),tmp1(idx),4),time_obs); % do polyfit
  idx = ~isnan(fit);
  obs_ts_em_hs        = fit-nanmean(fit(refstart-time_obs(1)+1:refende-time_obs(1)+1));  % do anomalies
  obs_residual        = am(obs_raw)-fit;

  obs_raw = squeeze(ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1.185001-202012_ts.nc'],'temperature'));
  obs = rm(am(obs_raw),wl);
  obs = obs-nanmean(obs(refstart-1850+1:refende-1850+1));

  % -- load GMST (need always)
  time_tas_obs = 1850:2020;
  tas_obs_raw = squeeze(ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1.185001-202012_ts.nc'],'temperature'));
  tas_obs = rm(am(tas_obs_raw),wl);
  ref_pi = nanmean(tas_obs(1850-1850+1:1899-1850+1)); % preindustrial reference
  ref = nanmean(tas_obs(refstart-1850+1:refende-1850+1));
  tas_obs_wtd = ref-ref_pi; % wtd = warming-to-date
  tas_obs = tas_obs-nanmean(tas_obs(refstart-1850+1:refende-1850+1));
end

if refende <= 2020
  % -- observations
  % -- spatial average of observations (needed for anything non-global mean)
  if strcmp(vari,'tas') == 1
    clear('tmp')
    obs_name = 'BEST';
    time_obs = 1850:2022;
    tmp0 = ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1_g025.nc'],'temperature');
    for i = 1:length(tmp0)
      tmp(i)     = nansum(nansum(squeeze(tmp0(jj,ii,i)).*weights));
    end
    if strcmp(seas,'annual')==1
      obs    = am(tmp);
    else
      obs    = seasmean(tmp,seas);
    end
    obs = obs-nanmean(obs(refstart-time_obs(1)+1:refende-time_obs(1)+1));
  else
    clear('tmp')
    obs_name = 'GPCC';
    time_obs = 1891:2016;
    tmp0 = ncread([pathin 'observations/precipitation/gpcc/precip.mon.mean.189101-201612_g025.nc'],'precip');
    for i = 1:length(tmp0)
      tmp(i)     = nansum(nansum(squeeze(tmp0(jj,ii,i)).*weights));
    end
    if strcmp(seas,'annual')==1
      obs    = am(tmp);
    else
      obs    = seasmean(tmp,seas);
    end
    obs = obs-nanmean(obs(refstart-time_obs(1)+1:refende-time_obs(1)+1));
  end
  obs_f = polyval(polyfit(time_obs,obs,4),time_obs); % polynomial fit
  obs_f = obs_f-nanmean(obs_f(refstart-time_obs(1)+1:refende-time_obs(1)+1));
  obs_r = obs-obs_f;
  obs_var = nanvar(rm(obs_r,wl))
end

% -- load CMIP6 data --
cmip6_ts_em  = NaN(length(models_cmip6),length(time));
scen        = scen_cmip6;
pathin_tmp  = [pathin 'cmip6-ng/' vari '/'];
for sc = 1:length(scen)
  fid = fopen(['~/Dropbox/tmp.txt'], 'w');
  clear('ne')
  for m = 1:length(models_cmip6)
    ['scenario = ' scen{sc} ' / model = ' models_cmip6{m} ]
    files   = dir([pathin 'cmip6-ng/' vari '/' vari '_mon_' models_cmip6{m} '_' scen{sc} '_r*_g025.nc']);
    ne(m) = 1;
    for e = 1:ne(m) %1:length(ensmem_cmip6(m))
      tmp0    = ncread([pathin_tmp files(e).name],[vari])*f; % yyy
      fprintf(fid, '%s\n', files(e).name);
      for i = 1:length(tmp0)
        if strcmp(region,'AA')==1
          tmp1(i)       = nansum(nansum(squeeze(tmp0(jj1,ii1,i)).*weights1)); % do area-weighted mean
          tmp2(i)       = nansum(nansum(squeeze(tmp0(jj2,ii2,i)).*weights2)); % do area-weighted mean
          tmp(i)     = tmp1(i)-tmp2(i);
        else
          tmp(i)     = nansum(nansum(squeeze(tmp0(jj,ii,i)).*weights));
        end
      end
      tmp     = tmp((start-start_cmip6)*12+1:(ende-start_cmip6+1)*12);
      if strcmp(seas,'annual')==1
        tmp    = am(tmp);
      else
        tmp    = seasmean(tmp,seas);
      end
      eval(['cmip6_' scen{sc} '_ts_raw(m,e,:) = tmp;'])
      ref = nanmean(tmp(refstart-start+1:refende-start+1));
      if strcmp(vari,'pr')==1
        eval(['cmip6_' scen{sc} '_ts(m,:)   = ((tmp-ref)/ref)*100;'])
      else
        eval(['cmip6_' scen{sc} '_ts(m,:)   = tmp-ref;'])
      end
    end
    % -- H&S style calculation --
    eval(['tmp1 = squeeze(cmip6_' scen{sc} '_ts_raw(m,1,:));'])
    idx = ~isnan(tmp1);
    fit = NaN(size(tmp1));
    fit(idx) = polyval(polyfit(time(idx),tmp1(idx)',4),time(idx)); % polynomial fit
    idx = ~isnan(fit);
    ref = nanmean(fit(refstart-start+1:refende-start+1));
    if strcmp(vari,'pr')==1
      eval(['cmip6_' scen{sc} '_ts_em(m,:)   = ((fit-ref)/ref)*100;'])
    else
      eval(['cmip6_' scen{sc} '_ts_em(m,:)   = fit-ref;'])
    end
    tmp1 = tmp1';
    residual  = tmp1-fit' + nanmean(fit);
    if strcmp(vari,'pr')==1
      ref       = nanmean(residual);
      residual  = (residual/ref)*100;
    end
    tmp2 = rm(residual,wl);
    eval(['cmip6_' scen{sc} '_ts_residual(m,:)  = tmp2;'])
    eval(['cmip6_' scen{sc} '_ts_var1(m)        = nanvar(tmp2);'])
    eval(['cmip6_' scen{sc} '_ts_var2(m)        = nanvar(tmp2(1:refende-start+1));'])
    eval(['cmip6_' scen{sc} '_ts_var3(m)        = nanvar(tmp2(refende-start+1+1:end));'])
  end
  fclose(fid);
end


% -- 2nd variable --
if refende <= 2020
  % -- Observations --
  % -- spatial average of observations (needed for anything non-global mean)
  if strcmp(vari2,'siconc')==0
    if strcmp(vari,'tas') == 1
      clear('tmp')
      obs2_name = 'BEST';
      time_obs2 = 1850:2022;
      tmp0 = squeeze(ncread([pathin 'observations/temperature/best/Land_and_Ocean_LatLong1.185001-202212_g025.nc'],'temperature'));
      for i = 1:length(tmp0)
        tmp(i)     = nansum(nansum(squeeze(tmp0(jj,ii,i)).*weights));
      end
      if strcmp(seas2,'annual')==1
        obs2    = am(tmp);
      else
        obs2    = seasmean(tmp,seas2);
      end
      obs2 = obs2-nanmean(obs2(refstart-time_obs2(1)+1:refende-time_obs2(1)+1));
    else
      clear('tmp')
      obs_name = 'GPCC';
      time_obs = 1891:2016;
      tmp0 = ncread([pathin 'observations/precipitation/gpcc/precip.mon.mean.189101-201612_g025.nc'],'precip');
      for i = 1:length(tmp0)
        tmp(i)     = nansum(nansum(squeeze(tmp0(jj,ii,i)).*weights));
      end
      if strcmp(seas2,'annual')==1
        obs2    = am(tmp);
      else
        obs2    = seasmean(tmp,seas2);
      end
      obs2 = obs2-nanmean(obs(refstart-1891+1:refende-1891+1));
    end
    obs2_f = polyval(polyfit(time_obs2,obs2,4),time_obs2); % polynomial fit
    obs2_f = obs2_f-nanmean(obs2_f(refstart-time_obs2(1)+1:refende-time_obs2(1)+1));
    obs2_r = obs2-obs2_f;
    obs2_var = nanvar(rm(obs2_r,wl));
  end
end

% -- CMIP6 models --
scen        = scen_cmip6;
pathin_tmp  = [pathin 'cmip6-ng/' vari2 '/'];
for sc = 1:length(scen)
  fid = fopen(['~/Dropbox/tmp.txt'], 'w');
  clear('ne')
  for m = 1:length(models_cmip6)
    ['scenario = ' scen{sc} ' / model = ' models_cmip6{m} ]
    if strcmp(vari2,'siconc')==1
      files   = dir([pathin 'cmip6-ng/' vari2 '/' vari2 '_mon_' models_cmip6{m} '_' scen{sc} '_r*_g025_NHts.nc']);
    else
      files   = dir([pathin 'cmip6-ng/' vari2 '/' vari2 '_mon_' models_cmip6{m} '_' scen{sc} '_r*_g025.nc']);
    end
    % if strcmp(scen{sc},'ssp585')==1
    %   ne(m) = length(files);
    % else
      ne(m) = 1;
    % end
    for e = 1:ne(m) %1:length(ensmem_cmip6(m))
      tmp0    = ncread([pathin_tmp files(e).name],[vari2])*f; % yyy
      fprintf(fid, '%s\n', files(e).name);
      if strcmp(vari2,'siconc')==1
        tmp = tmp0((start-start_cmip6)*12+1:(ende-start_cmip6+1)*12) * 1e-14;
        if strcmp(seas2,'S')==1
          eval(['cmip6_' scen{sc} '_ts_vari2(m,:)   = tmp(9:12:end);'])
        end
      else
        for i = 1:length(tmp0)
          if lnd_only == 1
            tmp(i)     = nansum(nansum(squeeze(tmp0(jj,ii,i)).*weights_lnd));
          else
            tmp(i)     = nansum(nansum(squeeze(tmp0(jj,ii,i)).*weights));
          end
        end
        tmp     = tmp((start-start_cmip6)*12+1:(ende-start_cmip6+1)*12);
        if strcmp(seas2,'annual')==1
          tmp    = am(tmp);
        else
          tmp    = seasmean(tmp,seas2);
        end
        eval(['cmip6_' scen{sc} '_ts_vari2_raw(m,e,:) = tmp;'])
        ref = nanmean(tmp(refstart-start+1:refende-start+1));
        if strcmp(vari2,'pr')==1
          eval(['cmip6_' scen{sc} '_ts_vari2(m,:)   = ((tmp-ref)/ref)*100;'])
        else
          eval(['cmip6_' scen{sc} '_ts_vari2(m,:)   = tmp-ref;'])
        end
      end
    end
  end
  fclose(fid);
end


% -- load CMIP6 piControl data --
clear('tmp','tmp1','tmp2')
cmip6_piControl_ts_raw  = NaN(length(models_cmip6),252);
cmip6_piControl_ts      = NaN(length(models_cmip6),252);
for m = 1:length(models_cmip6)
  ['model = ' models_cmip6{m}]
  filein  = dir([pathin 'cmip6-ng/' vari '/' vari '_mon_' models_cmip6{m} '_piControl_*_g025.nc']);
  tmp0    = ncread([filein(1).folder '/' filein(1).name],[vari])*f;
  tmp0    = tmp0(:,:,end-252*12+1:end);
  for i = 1:length(tmp0)
    if strcmp(region,'AA')==1
      tmp1(i)    = nansum(nansum(squeeze(tmp0(jj1,ii1,i)).*weights1)); % do area-weighted mean
      tmp2(i)    = nansum(nansum(squeeze(tmp0(jj2,ii2,i)).*weights2)); % do area-weighted mean
      tmp(i)     = tmp1(i)-tmp2(i);
    else
      tmp(i)     = nansum(nansum(squeeze(tmp0(jj,ii,i)).*weights));
    end
  end
  if strcmp(seas,'annual')==1
    tmp    = am(tmp);
  else
    tmp    = seasmean(tmp,seas);
  end
  cmip6_piControl_ts_raw(m,:) = tmp;
  if wl > 1
    tmp = rm(tmp,wl);
  else
    tmp = tmp;
  end
  ref = nanmean(tmp);
  if strcmp(vari,'pr')==1
    cmip6_piControl_ts(m,:)   = ((tmp-ref)/ref)*100;
  else
    cmip6_piControl_ts(m,:)   = tmp-ref;
  end
  idx = ~isnan(cmip6_piControl_ts(m,:));
  cmip6_piControl_ts_var(m)       = nanvar(detrend(cmip6_piControl_ts(m,idx)));
  cmip6_piControl_ts_vari_std(m)  = nanstd(cmip6_piControl_ts(m,idx));
end
if strcmp(vari2,'siconc')==0
  % -- vari2
  cmip6_piControl_ts_raw_vari2  = NaN(length(models_cmip6),252);
  cmip6_piControl_ts_vari2      = NaN(length(models_cmip6),252);
  for m = 1:length(models_cmip6)
    ['model = ' models_cmip6{m}]
    filein  = dir([pathin 'cmip6-ng/' vari2 '/' vari2 '_mon_' models_cmip6{m} '_piControl_*_g025.nc']);
    tmp0    = ncread([filein(1).folder '/' filein(1).name],[vari2])*f;
    tmp0    = tmp0(:,:,end-252*12+1:end);
    for i = 1:length(tmp0)
      if strcmp(region,'AA')==1
        tmp1(i)    = nansum(nansum(squeeze(tmp0(jj1,ii1,i)).*weights1)); % do area-weighted mean
        tmp2(i)    = nansum(nansum(squeeze(tmp0(jj2,ii2,i)).*weights2)); % do area-weighted mean
        tmp(i)     = tmp1(i)-tmp2(i);
      else
        tmp(i)     = nansum(nansum(squeeze(tmp0(jj,ii,i)).*weights));
      end
    end
    if strcmp(seas2,'annual')==1
      tmp    = am(tmp);
    else
      tmp    = seasmean(tmp,seas2);
    end
    cmip6_piControl_ts_raw_vari2(m,:) = tmp;
    if wl > 1
      tmp = rm(tmp,wl);
    else
      tmp = tmp;
    end
    ref = nanmean(tmp);
    if strcmp(vari,'pr')==1
      cmip6_piControl_ts_vari2(m,:)   = ((tmp-ref)/ref)*100;
    else
      cmip6_piControl_ts_vari2(m,:)   = tmp-ref;
    end
    idx = ~isnan(cmip6_piControl_ts_vari2(m,:));
    cmip6_piControl_ts_vari2_var(m)   = nanvar(detrend(cmip6_piControl_ts_vari2(m,idx)));
  end
  for m = 1:length(models_cmip6)
    cmip6_piControl_ts_vari_vari2_cor(m) = nanmean(corr(cmip6_piControl_ts(m,idx)',cmip6_piControl_ts_vari2(m,idx)'));
    cmip6_piControl_ts_vari2_std(m)      = nanstd(cmip6_piControl_ts_vari2(m,idx));
  end
end


%% -- SCALING CALCULATIONS -----------
% -- CMIP6 scen uncertainty with 4 scenarios, but scaled ("squeezed") to reflect the fact that some scenarios have become less likely in the recent decade:
% -- motivated by Moore et al., 2022, Nature
cmip6_scen = [...
nanmean(cmip6_ssp585_ts_em,1);...
nanmean(cmip6_ssp370_ts_em,1);...
nanmean(cmip6_ssp245_ts_em,1);...
nanmean(cmip6_ssp126_ts_em,1);...
nanmean(cmip6_ssp119_ts_em,1)];
% --
n = length(time);
% -- weights
f1 = .60; % ssp585
f2 = .65; % ssp370
f3 = .8; % ssp245
f4 = 1.1; % ssp126
f5 = 1; % ssp119
sc1 = time*0;
sc1(1:refende-start) = 1;
sc2 = sc1;
sc3 = sc1;
sc4 = sc1;
sc5 = sc1;
sc1(refende-start+1:n) = linspace(1,f1,length(refende-start+1:n));
sc2(refende-start+1:n) = linspace(1,f2,length(refende-start+1:n));
sc3(refende-start+1:n) = linspace(1,f3,length(refende-start+1:n));
sc4(refende-start+1:n) = linspace(1,f4,length(refende-start+1:n));
sc5(refende-start+1:n) = linspace(1,f5,length(refende-start+1:n));
cmip6_scen_moore = [...
nanmean(cmip6_ssp585_ts_em,1).*sc1;...
nanmean(cmip6_ssp370_ts_em,1).*sc2;...
nanmean(cmip6_ssp245_ts_em,1).*sc3;...
nanmean(cmip6_ssp126_ts_em,1).*sc4;...
nanmean(cmip6_ssp119_ts_em,1).*sc5];
for i = 1:length(scen_cmip6)
  cmip6_scen_moore(i,1:refstart-start) = cmip6_scen_moore(3,1:refstart-start);
end

% -- observed warming since preindustrial (today = reference period)
start_pi = 1850;
ende_pi  = 1900;
if refende <= 2020 && strcmp(vari,'tas')==1
  if strcmp(region,'global')==1
    offset = nanmean(obs_raw_am(refstart-time_obs(1)+1:refende-time_obs(1))) - nanmean(obs_raw_am(start_pi-time_obs(1)+1:ende_pi-time_obs(1)))
  else
    offset = nanmean(obs(refstart-time_obs(1)+1:refende-time_obs(1))) - nanmean(obs(start_pi-time_obs(1)+1:ende_pi-time_obs(1)))
  end
else
  offset = 0;
end


% -- Qasmi+Ribes (2022) local constraint --
pathin_tmp  = '/Users/fl439/Dropbox/work/cmip6_ipcc_ar6_wg1_assessed_ranges/';
clear('tmp1','tmp2','ssp585','ssp585_loc','ssp585_loc_ts')
tmp1(:,:,:,1)     = ncread([pathin_tmp 'K4CMIP_fullfield_tblend_CMIP6_histssp585_ann_g025_ok.nc'],'q05_uncons_loc');
tmp1(:,:,:,2)     = ncread([pathin_tmp 'K4CMIP_fullfield_tblend_CMIP6_histssp585_ann_g025_ok.nc'],'avg_uncons_loc');
tmp1(:,:,:,3)     = ncread([pathin_tmp 'K4CMIP_fullfield_tblend_CMIP6_histssp585_ann_g025_ok.nc'],'q95_uncons_loc');
tmp1_loc(:,:,:,1) = ncread([pathin_tmp 'K4CMIP_fullfield_tblend_CMIP6_histssp585_ann_g025_ok.nc'],'q05_cons_loc_glo');
tmp1_loc(:,:,:,2) = ncread([pathin_tmp 'K4CMIP_fullfield_tblend_CMIP6_histssp585_ann_g025_ok.nc'],'avg_cons_loc_glo');
tmp1_loc(:,:,:,3) = ncread([pathin_tmp 'K4CMIP_fullfield_tblend_CMIP6_histssp585_ann_g025_ok.nc'],'q95_cons_loc_glo');
% -- need to remap from -180:180 to 0:360
ssp585      = tmp1;
ssp585_loc  = tmp1_loc;
ssp585(1:72,:,:,:)        = tmp1(73:end,:,:,:);
ssp585(73:end,:,:,:)      = tmp1(1:72,:,:,:);
ssp585_loc(1:72,:,:,:)    = tmp1_loc(73:end,:,:,:);
ssp585_loc(73:end,:,:,:)  = tmp1_loc(1:72,:,:,:);
j = 1;
for i = 101:length(ssp585_loc)-1
  ssp585_ts(j,1)      = nansum(nansum(squeeze(ssp585(jj,ii,i,1)).*weights))';
  ssp585_ts(j,2)      = nansum(nansum(squeeze(ssp585(jj,ii,i,2)).*weights))';
  ssp585_ts(j,3)      = nansum(nansum(squeeze(ssp585(jj,ii,i,3)).*weights))';
  ssp585_loc_ts(j,1)  = nansum(nansum(squeeze(ssp585_loc(jj,ii,i,1)).*weights))';
  ssp585_loc_ts(j,2)  = nansum(nansum(squeeze(ssp585_loc(jj,ii,i,2)).*weights))';
  ssp585_loc_ts(j,3)  = nansum(nansum(squeeze(ssp585_loc(jj,ii,i,3)).*weights))';
  j = j+1;
end
% -- same reference period as cmip6 simulations:
ssp585_ts     = ssp585_ts-nanmean(ssp585_ts(refstart-start+1:refende-start,2));
ssp585_loc_ts = ssp585_loc_ts-nanmean(ssp585_loc_ts(refstart-start+1:refende-start,2));
% -- how much is uncertainty reduced in Qasmi&Ribes?
% -- absolute warming (ssp585):
ssp585_loc_ts(101,2)./ssp585_ts(101,2) % at 2050
ssp585_loc_ts(end,2)./ssp585_ts(end,2) % at 2099
% -- range of warming (ssp585):
(ssp585_loc_ts(101,3)-ssp585_loc_ts(101,1))./(ssp585_ts(101,3)-ssp585_ts(101,1)) % at 2050
(ssp585_loc_ts(end,3)-ssp585_loc_ts(end,1))./(ssp585_ts(end,3)-ssp585_ts(end,1)) % at 2099
% -- determine scaling from ssp585, to be used for other scenarios
ratio_offset  = 300;
ssp585_loc_ts_mean_ratio = (ssp585_loc_ts(:,2)+ratio_offset)./(ssp585_ts(:,2)+ratio_offset);
range_p95_mean      = ssp585_ts(:,3)-ssp585_ts(:,2);
range_p95_mean_loc  = ssp585_loc_ts(:,3)-ssp585_loc_ts(:,2);
range_p5_mean       = ssp585_ts(:,2)-ssp585_ts(:,1);
range_p5_mean_loc   = ssp585_loc_ts(:,2)-ssp585_loc_ts(:,1);
ssp585_loc_ts_p95_mean_ratio  = range_p95_mean_loc./range_p95_mean;
ssp585_loc_ts_p5_mean_ratio   = range_p5_mean_loc./range_p5_mean;
% -- use the Qasmi&Rhibes ssp585 constraint as scaling for other scenarios --
ssp585_loc2_ts(:,1) = (nanmean(cmip6_ssp585_ts_em)+ratio_offset).*ssp585_loc_ts_mean_ratio'-ratio_offset - ((nanmean(cmip6_ssp585_ts_em)-prctile(cmip6_ssp585_ts_em,5)).*ssp585_loc_ts_p5_mean_ratio');
ssp585_loc2_ts(:,2) = (nanmean(cmip6_ssp585_ts_em)+ratio_offset).*ssp585_loc_ts_mean_ratio'-ratio_offset;
ssp585_loc2_ts(:,3) = (prctile(cmip6_ssp585_ts_em,95)-nanmean(cmip6_ssp585_ts_em)).*ssp585_loc_ts_p95_mean_ratio' + (nanmean(cmip6_ssp585_ts_em)+ratio_offset).*ssp585_loc_ts_mean_ratio'-ratio_offset;
% --
ssp370_loc_ts_mean_ratio = ssp585_loc_ts_mean_ratio;
ssp370_loc_ts_mean_ratio(1:refstart-start+10) = 1;
ssp370_loc_ts_mean_ratio(refstart-start+10+1:end) = linspace(ssp585_loc_ts_mean_ratio(refstart-start+10+1),ssp585_loc_ts_mean_ratio(end)*1.0004,2099-refstart-10+1);
ssp370_loc2_ts(:,1) = (nanmean(cmip6_ssp370_ts_em)+ratio_offset).*ssp370_loc_ts_mean_ratio'-ratio_offset - ((nanmean(cmip6_ssp370_ts_em)-prctile(cmip6_ssp370_ts_em,5)).*ssp585_loc_ts_p5_mean_ratio');
ssp370_loc2_ts(:,2) = (nanmean(cmip6_ssp370_ts_em)+ratio_offset).*ssp370_loc_ts_mean_ratio'-ratio_offset;
ssp370_loc2_ts(:,3) = (prctile(cmip6_ssp370_ts_em,95)-nanmean(cmip6_ssp370_ts_em)).*ssp585_loc_ts_p95_mean_ratio' + (nanmean(cmip6_ssp370_ts_em)+ratio_offset).*ssp370_loc_ts_mean_ratio'-ratio_offset;
% --
ssp245_loc_ts_mean_ratio = ssp585_loc_ts_mean_ratio;
ssp245_loc_ts_mean_ratio(1:refstart-start+10) = 1;
ssp245_loc_ts_mean_ratio(refstart-start+10+1:end) = linspace(ssp585_loc_ts_mean_ratio(refstart-start+10+1),ssp585_loc_ts_mean_ratio(end)*1.0008,2099-refstart-10+1);
ssp245_loc2_ts(:,1) = (nanmean(cmip6_ssp245_ts_em)+ratio_offset).*ssp245_loc_ts_mean_ratio'-ratio_offset - ((nanmean(cmip6_ssp245_ts_em)-prctile(cmip6_ssp245_ts_em,5)).*ssp585_loc_ts_p5_mean_ratio');
ssp245_loc2_ts(:,2) = (nanmean(cmip6_ssp245_ts_em)+ratio_offset).*ssp245_loc_ts_mean_ratio'-ratio_offset;
ssp245_loc2_ts(:,3) = (prctile(cmip6_ssp245_ts_em,95)-nanmean(cmip6_ssp245_ts_em)).*ssp585_loc_ts_p95_mean_ratio' + (nanmean(cmip6_ssp245_ts_em)+ratio_offset).*ssp245_loc_ts_mean_ratio'-ratio_offset;
% --
ssp126_loc_ts_mean_ratio = ssp585_loc_ts_mean_ratio;
ssp126_loc_ts_mean_ratio(1:refstart-start+10) = 1;
ssp126_loc_ts_mean_ratio(refstart-start+10+1:end) = linspace(ssp585_loc_ts_mean_ratio(refstart-start+10+1),ssp585_loc_ts_mean_ratio(end)*1.0012,2099-refstart-10+1);
ssp126_loc2_ts(:,1) = (nanmean(cmip6_ssp126_ts_em)+ratio_offset).*ssp126_loc_ts_mean_ratio'-ratio_offset - ((nanmean(cmip6_ssp126_ts_em)-prctile(cmip6_ssp126_ts_em,5)).*ssp585_loc_ts_p5_mean_ratio');
ssp126_loc2_ts(:,2) = (nanmean(cmip6_ssp126_ts_em)+ratio_offset).*ssp126_loc_ts_mean_ratio'-ratio_offset;
ssp126_loc2_ts(:,3) = (prctile(cmip6_ssp126_ts_em,95)-nanmean(cmip6_ssp126_ts_em)).*ssp585_loc_ts_p95_mean_ratio' + (nanmean(cmip6_ssp126_ts_em)+ratio_offset).*ssp126_loc_ts_mean_ratio'-ratio_offset;
% --
if length(scen_cmip6)==5
  ssp119_loc_ts_mean_ratio = ssp585_loc_ts_mean_ratio;
  ssp119_loc_ts_mean_ratio(1:refstart-start+10) = 1;
  ssp119_loc_ts_mean_ratio(refstart-start+10+1:end) = linspace(ssp585_loc_ts_mean_ratio(refstart-start+10+1),ssp585_loc_ts_mean_ratio(end)*1.0016,2099-refstart-10+1);
  ssp119_loc2_ts(:,1) = (nanmean(cmip6_ssp119_ts_em)+ratio_offset).*ssp119_loc_ts_mean_ratio'-ratio_offset - ((nanmean(cmip6_ssp119_ts_em)-prctile(cmip6_ssp119_ts_em,5)).*ssp585_loc_ts_p5_mean_ratio');
  ssp119_loc2_ts(:,2) = (nanmean(cmip6_ssp119_ts_em)+ratio_offset).*ssp119_loc_ts_mean_ratio'-ratio_offset;
  ssp119_loc2_ts(:,3) = (prctile(cmip6_ssp119_ts_em,95)-nanmean(cmip6_ssp119_ts_em)).*ssp585_loc_ts_p95_mean_ratio' + (nanmean(cmip6_ssp119_ts_em)+ratio_offset).*ssp119_loc_ts_mean_ratio'-ratio_offset;
end
% -- scale with new scenario uncertainty:
ssp126_loc2_ts_moore = ssp126_loc2_ts.*sc4';
ssp245_loc2_ts_moore = ssp245_loc2_ts.*sc3';
ssp370_loc2_ts_moore = ssp370_loc2_ts.*sc2';
ssp585_loc2_ts_moore = ssp585_loc2_ts.*sc1';
if length(scen_cmip6)==5
  ssp119_loc2_ts_moore = ssp119_loc2_ts.*sc5';
end


% -- UNCERTAINTIES --
% -- CMIP6 scen uncertainty with 5 scenarios:
scen_u_cmip6      = nanvar([nanmean(cmip6_ssp585_ts_em,1); nanmean(cmip6_ssp370_ts_em,1); nanmean(cmip6_ssp245_ts_em,1); nanmean(cmip6_ssp126_ts_em,1); nanmean(cmip6_ssp119_ts_em,1)]); % using 5 scenarios
% -- CMIP6 scen uncertainty after scaling by Moore 2022:
scen_u_cmip6_moore = nanvar(cmip6_scen_moore);
% -- CMIP6 scen uncertainty after scaling model uncertainty with Qasmi:
scen_u_cmip6_loc2 = nanvar([ssp119_loc2_ts(:,2),ssp126_loc2_ts(:,2),ssp245_loc2_ts(:,2),ssp370_loc2_ts(:,2),ssp585_loc2_ts(:,2)]');
% -- CMIP6 scen uncertainty after scaling model uncertainty with Qasmi and scenario uncertainty with Moore:
scen_u_cmip6_loc2_moore = nanvar([ssp119_loc2_ts_moore(:,2)'; ssp126_loc2_ts_moore(:,2)'; ssp245_loc2_ts_moore(:,2)'; ssp370_loc2_ts_moore(:,2)'; ssp585_loc2_ts_moore(:,2)']);

model_u_cmip6       = nanvar(cmip6_ssp245_ts_em,1);
model_u_cmip6_moore = nanvar(cmip6_ssp245_ts_em.*sc3,1);

% -- HS09: calculate the mean across model uncertainties from the different scenarios):
model_u_cmip6_hs        = nanmean([nanvar(cmip6_ssp585_ts_em,1); nanvar(cmip6_ssp370_ts_em,1); nanvar(cmip6_ssp245_ts_em,1); nanvar(cmip6_ssp126_ts_em,1)]);
model_u_cmip6_hs_3scen  = nanmean([nanvar(cmip6_ssp370_ts_em,1); nanvar(cmip6_ssp245_ts_em,1); nanvar(cmip6_ssp126_ts_em,1)]);

% -- model uncertainty from constrained ranges of Qasmi+Ribes (2022)
% -- derive stddev from confidence interval following https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm#:~:text=The%20standard%20deviation%20for%20each,should%20be%20replaced%20by%205.15.
% -- actual constrained projections from paper
tmp = ssp585_loc_ts-nanmean(ssp585_loc_ts(refstart-1950+1:refende-1950,2)); % make sure multi-model mean is zero at reference period
model_u_cmip6_loc = (((tmp(:,3)-tmp(:,1))/3.29).^2)';
% -- ad hoc scaled ssp245 projections (use actual ssp245 once they send it)
tmp = ssp245_loc2_ts-nanmean(ssp245_loc2_ts(refstart-1950+1:refende-1950,2)); % make sure multi-model mean is zero at reference period
model_u_cmip6_loc2 = (((tmp(:,3)-tmp(:,1))/3.29).^2)';
% -- scaled with scenario uncertainty scaling from Moore
tmp = ssp245_loc2_ts_moore-nanmean(ssp245_loc2_ts_moore(refstart-1950+1:refende-1950,2)); % make sure multi-model mean is zero at reference period
model_u_cmip6_loc2_moore = (((tmp(:,3)-tmp(:,1))/3.29).^2)';

int_u_cmip6_mean  = nanmean(cmip6_ssp585_ts_var1);
int_u_cmip6_max   = max(cmip6_ssp585_ts_var1);
int_u_cmip6_min   = min(cmip6_ssp585_ts_var1);

total_u_cmip6_mean        = scen_u_cmip6 + model_u_cmip6 + int_u_cmip6_mean;
total_u_cmip6_mean_moore        = scen_u_cmip6_moore + model_u_cmip6_moore + int_u_cmip6_mean;
total_u_cmip6_mean_loc2         = scen_u_cmip6_loc2 + model_u_cmip6_loc2 + int_u_cmip6_mean;
total_u_cmip6_mean_loc2_moore   = scen_u_cmip6_loc2_moore + model_u_cmip6_loc2_moore + int_u_cmip6_mean;
total_u_cmip6_mean_hs     = scen_u_cmip6 + model_u_cmip6_hs + int_u_cmip6_mean;

% -- fractional uncertainties --
scen_u_frac_cmip6     = (scen_u_cmip6./total_u_cmip6_mean)*100;
scen_u_frac_cmip6_moore       = (scen_u_cmip6_moore./total_u_cmip6_mean_moore)*100;
scen_u_frac_cmip6_loc2        = (scen_u_cmip6_loc2./total_u_cmip6_mean_loc2)*100;
scen_u_frac_cmip6_loc2_moore  = (scen_u_cmip6_loc2_moore./total_u_cmip6_mean_loc2_moore)*100;
model_u_frac_cmip6    = (model_u_cmip6./total_u_cmip6_mean)*100;
model_u_frac_cmip6_moore        = (model_u_cmip6_moore./total_u_cmip6_mean_moore)*100;
model_u_frac_cmip6_loc2         = (model_u_cmip6_loc2./total_u_cmip6_mean_loc2)*100;
model_u_frac_cmip6_loc2_moore   = (model_u_cmip6_loc2_moore./total_u_cmip6_mean_loc2_moore)*100;
int_u_frac_cmip6_mean     = (int_u_cmip6_mean./total_u_cmip6_mean)*100;
int_u_frac_cmip6_mean_moore   = (int_u_cmip6_mean./total_u_cmip6_mean_moore)*100;
int_u_frac_cmip6_mean_loc2        = (int_u_cmip6_mean./total_u_cmip6_mean_loc2)*100;
int_u_frac_cmip6_mean_loc2_moore  = (int_u_cmip6_mean./total_u_cmip6_mean_loc2_moore)*100;



% --- PLOTTING ---------------------------------------------------------------
close all

cols = [255 0 0;
      255 160 16;
      255 224 32;
      0 192 0;
      80 208 255;
      0 32 255;
      160 32 255]/255;

hs09_cols = ...
[53 74 161;...
255 110 4;...
0 127 60]/255;
hs09_cols_light = ...
[164 180 245;...
252 210 179;...
172 232 200]/255;

rcp_cols = [217,37,42;...
         232,126,63;...
         133,177,212;...
         51,74,141]/255;
rcp_cols_light = ...
[235,144,115;...
243,183,136;...
189,211,227;...
138,141,185]/255;
rcp_cols_light2 = (1-rcp_cols)*.25+rcp_cols;

ssp_cols = [25,162,192;...
         42,55,89;...
         223,137,63;...
         201,54,56,;...
         131,38,40]/255;

if wl > 1
xlim  = [2015 2099-round(wl/2)];
xlim0 = [1950+round(wl/2) 2099-round(wl/2)];
else
xlim  = [refende+1 2099-round(wl/2)];
xlim0 = [1950+round(wl/2) 2099-round(wl/2)];
end
if strcmp(vari,'tas')==1
  yincr = 1;
  if refstart == 1995 || refstart == 2001
    ylim = [-1.25 5.5]*f;
  else
    ylim = [0 5.5];
  end
else % vari = pr
  yincr = 2;
  ylim = [-2.9 10];
  if strcmp(region,'Mediterranean')==1
    ylim = [-10 10];
  end
  if strcmp(region,'UCRB')==1
    ylim = [-10 20];
  end
end

tl = 0.015; % tick length
sh = 0.05; % horizontal space between panels
sv = 0.04; % vertical space between panels


% -- plot --------
close all

end_year = 2050; % 2050 2095
xlim = [2011 end_year];
if strcmp(region,'global')==1
if end_year == 2050
  ylim = [-.2 2.5];
else
  ylim = [-.2 7];
end
elseif strcmp(region,'arctic')==1 || strcmp(region,'Hudson_Bay')==1
if end_year == 2050
  ylim = [-.3*2.2 3.75*2.2];
else
  ylim = [-.6*2.2 8*2.2];
end
else
if end_year == 2050
  ylim = [-.3 3.75];
else
  ylim = [-.6 8];
end
end

figure1 = figure;
set(figure1, 'units', 'centimeters', 'pos', [10 10 30 12]) % final size

% -- CMIP6:
subplot(1,3,1)
hold on
title([region ' ' vari ' ' seas ' ' num2str(wl) '-yr means'],'Interpreter','none')%,'FontSize',10)
% -- scenario scaling and model u scaling
jbfill(time,prctile(cmip6_ssp585_ts_em,95),prctile(cmip6_ssp585_ts_em,5),rcp_cols_light(1,:),'none',1,.2) % raw
jbfill(time,ssp585_loc2_ts_moore(:,3)',ssp585_loc2_ts_moore(:,1)',rcp_cols_light(1,:),'none',1,.6) % scen+model constrained
jbfill(time,prctile(cmip6_ssp119_ts_em,95),prctile(cmip6_ssp119_ts_em,5),rcp_cols_light(4,:),'none',1,.2) % raw
jbfill(time,ssp119_loc2_ts_moore(:,3)',ssp119_loc2_ts_moore(:,1)',rcp_cols_light(4,:),'none',1,.6) % scen+model constrained
% --
% -- plot respective multi-model means
hold on
h1 = plot(time,nanmean(cmip6_ssp585_ts_em),'Color',rcp_cols(1,:),'LineWidth',2);
plot(time,nanmean(cmip6_ssp119_ts_em),'Color',rcp_cols(4,:),'LineWidth',2)
plot(time,ssp119_loc2_ts(:,2),'--','Color',rcp_cols(4,:),'LineWidth',2)
plot(time,nanmean(cmip6_ssp119_ts_em.*sc5),':','Color',rcp_cols(4,:),'LineWidth',2)
h4 = plot(time,ssp585_loc2_ts_moore(:,2),'Color',rcp_cols(1,:),'LineStyle','-.','LineWidth',2);
plot(time,ssp119_loc2_ts_moore(:,2),'Color',rcp_cols(4,:),'LineStyle','-.','LineWidth',2)
if refende <= 2020
  plot(time_obs,obs,'Color',[.3 .3 .3],'LineWidth',1)
  plot(time_obs,rm(obs,10),'k','LineWidth',3)
end
set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'XLim',xlim,'TickLength',[tl tl],'Layer','top')
hline(0,'k')
box on
legend([h1(1) h4(1)],'Unconstrained',['Scenario and response uncertainty' char(10) 'constrained'],'Location','NorthWest','Interpreter','none')
legend boxoff
xlabel('Time (Year)')
ylabel('Temperature change (\circC)')
yyaxis right
set(gca,'YLim',ylim+offset,'YTick',[ceil(ylim(1)):yincr:10],'XLim',xlim)
ylabel(['relative to preindustrial (\circC)'])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';


% ------------
subplot(1,3,2)
hold on
title(['Temperature range at year ' num2str(end_year)])
['Temperature range at year ' num2str(end_year) ':']
os = 0; % .1-.05
j = end_year-start+1; % point in time
x  = 1;
if length(scen_cmip6)==5
  plot([x-os x-os],[prctile(cmip6_ssp119_ts_em(:,j),5) prctile(cmip6_ssp585_ts_em(:,j),95)],'k','LineWidth',10)
  ['Unconstrained = ' num2str([prctile(cmip6_ssp119_ts_em(:,j),5) prctile(cmip6_ssp585_ts_em(:,j),95)] + offset) ]
  range_unconstrained = prctile(cmip6_ssp585_ts_em(:,j),95) - prctile(cmip6_ssp119_ts_em(:,j),5);
else
  plot([x-os x-os],[prctile(cmip6_ssp126_ts_em(:,j),5) prctile(cmip6_ssp585_ts_em(:,j),95)],'k','LineWidth',10)
  ['Unconstrained = ' num2str([prctile(cmip6_ssp126_ts_em(:,j),5) prctile(cmip6_ssp585_ts_em(:,j),95)] + offset) ]
  range_unconstrained = prctile(cmip6_ssp585_ts_em(:,j),95) - prctile(cmip6_ssp126_ts_em(:,j),5);
end
x  = 1.5;
tmp1 = cmip6_ssp119_ts_em.*sc5;
tmp2 = cmip6_ssp585_ts_em.*sc1;
plot([x-os x-os],[prctile(tmp1(:,j),5) prctile(tmp2(:,j),95)],'k','LineWidth',10)
['Scenario uncertainty constrained = ' num2str([prctile(tmp1(:,j),5) prctile(tmp2(:,j),95)] + offset) ]
range_scen_constrained = prctile(tmp2(:,j),95) - prctile(tmp1(:,j),5);
x  = 2;
plot([x-os x-os],[ssp119_loc2_ts(j,1) ssp585_loc2_ts(j,3)],'k','LineWidth',10)
['Response uncertainty constrained = ' num2str([ssp119_loc2_ts(j,1) ssp585_loc2_ts(j,3)] + offset) ]
range_model_constrained = ssp585_loc2_ts(j,3) - ssp119_loc2_ts(j,1);
x  = 2.5;
p0 = plot([x-os x-os],[ssp119_loc2_ts_moore(j,1) ssp585_loc2_ts_moore(j,3)],'k','LineWidth',10);
['Scenario and Response uncertainty constrained = ' num2str([ssp119_loc2_ts_moore(j,1) ssp585_loc2_ts_moore(j,3)] + offset) ]
range_scen_model_constrained = ssp585_loc2_ts_moore(j,3) - ssp119_loc2_ts_moore(j,1);
% -
plot([2.5-os xlim(2)],[ssp585_loc2_ts_moore(j,3) ssp585_loc2_ts_moore(j,3)],'k:','LineWidth',1)
plot([2.5-os xlim(2)],[ssp119_loc2_ts_moore(j,1) ssp119_loc2_ts_moore(j,1)],'k:','LineWidth',1)
if strcmp(region,'global')==1
  text(.88,.9,'5-95% range of all projections','rotation',90)%,'FontWeight','bold')
elseif strcmp(region,'Central_Europe')==1
  text(.87,1.05,'5-95% range of all projections','rotation',90)%,'FontWeight','bold')
else
  text(.88,.9,'5-95% range of all projections','rotation',90)%,'FontWeight','bold')
end
set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],...
'XTick',[1:.5:2.5],'xticklabel',{'Unconstrained';'Scenario uncertainty constrained';'Response uncertainty constrained';'Scenario and response uncertainty constrained'},...
'TickLength',[tl tl],'Layer','top')
xlim = [.6 2.9];
set(gca,'XLim',xlim)
hline(0,'k')
box on
% --
yyaxis right
set(gca,'YLim',ylim+offset,'YTick',[ceil(ylim(1)):yincr:10],'XLim',xlim)
ylabel(['relative to preindustrial (\circC)'])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
% -- pie charts section
he = .12; % y-axis position
si = .14; % pie size?
axes('Position',[.378 he si si])
{'Response U', 'Int. var.', 'Scenario U'}
tmp = [model_u_frac_cmip6(j) int_u_frac_cmip6_mean(j) scen_u_frac_cmip6(j)]
p = pie(tmp, repmat({''},size(tmp)));
colormap(hs09_cols)

axes('Position',[.425 (he*(1+0.5*(1-(range_scen_constrained/range_unconstrained)))) si si*(range_scen_constrained/range_unconstrained)]) % scale pie chart size by total uncertainty
tmp = [model_u_frac_cmip6_moore(j) int_u_frac_cmip6_mean_moore(j) scen_u_frac_cmip6_moore(j)]
pie(tmp, repmat({''},size(tmp)))

axes('Position',[.472 (he*(1+0.5*(1-(range_model_constrained/range_unconstrained)))) si si*(range_model_constrained/range_unconstrained)]) % scale pie chart size by total uncertainty
tmp = [model_u_frac_cmip6_loc2(j) int_u_frac_cmip6_mean_loc2(j) scen_u_frac_cmip6_loc2(j)]
pie(tmp, repmat({''},size(tmp)))

axes('Position',[.518 (he*(1+0.5*(1-(range_scen_model_constrained/range_unconstrained)))) si si*(range_scen_model_constrained/range_unconstrained)]) % scale pie chart size by total uncertainty
tmp = [model_u_frac_cmip6_loc2_moore(j) int_u_frac_cmip6_mean_loc2_moore(j) scen_u_frac_cmip6_loc2_moore(j)]
pie(tmp, repmat({''},size(tmp)))

l = legend([p(3) p(1) p(5)],'Int. variability','Response','Scenario');
pos = get(l,'Position');
posx = 0.5;
posy = 0.8;
set(l,'Position',[posx posy pos(3) pos(4)]);
% -- end pie chart section


% ------------
subplot(1,3,3)
hold on
if lnd_only == 1 && strcmp(region,'global')==1
  title(['2nd constraint: ' seas2 ' ' vari2 ' ' region ' lnd'])
else
  title(['2nd constraint: ' seas2 ' ' vari2 ' ' region])
end
tmp_raw       = [cmip6_ssp119_ts_em(:,j-5:j+4); cmip6_ssp126_ts_em(:,j-5:j+4); cmip6_ssp245_ts_em(:,j-5:j+4); cmip6_ssp370_ts_em(:,j-5:j+4); cmip6_ssp585_ts_em(:,j-5:j+4)];
tmp_vari2_raw = [cmip6_ssp119_ts_vari2(:,j-5:j+4); cmip6_ssp126_ts_vari2(:,j-5:j+4); cmip6_ssp245_ts_vari2(:,j-5:j+4); cmip6_ssp370_ts_vari2(:,j-5:j+4); cmip6_ssp585_ts_vari2(:,j-5:j+4)];
% -- clobber all scenarios together:
tmp           = [nanmean(cmip6_ssp119_ts_em(:,j-5:j+4),2); nanmean(cmip6_ssp126_ts_em(:,j-5:j+4),2); nanmean(cmip6_ssp245_ts_em(:,j-5:j+4),2); nanmean(cmip6_ssp370_ts_em(:,j-5:j+4),2); nanmean(cmip6_ssp585_ts_em(:,j-5:j+4),2)];
tmp_vari2     = [nanmean(cmip6_ssp119_ts_vari2(:,j-5:j+4),2); nanmean(cmip6_ssp126_ts_vari2(:,j-5:j+4),2); nanmean(cmip6_ssp245_ts_vari2(:,j-5:j+4),2); nanmean(cmip6_ssp370_ts_vari2(:,j-5:j+4),2); nanmean(cmip6_ssp585_ts_vari2(:,j-5:j+4),2)];
% -- keep scenario categories:
tmp_cat       = [nanmean(cmip6_ssp119_ts_em(:,j-5:j+4),2) nanmean(cmip6_ssp126_ts_em(:,j-5:j+4),2) nanmean(cmip6_ssp245_ts_em(:,j-5:j+4),2) nanmean(cmip6_ssp370_ts_em(:,j-5:j+4),2) nanmean(cmip6_ssp585_ts_em(:,j-5:j+4),2)];
tmp_vari2_cat = [nanmean(cmip6_ssp119_ts_vari2(:,j-5:j+4),2) nanmean(cmip6_ssp126_ts_vari2(:,j-5:j+4),2) nanmean(cmip6_ssp245_ts_vari2(:,j-5:j+4),2) nanmean(cmip6_ssp370_ts_vari2(:,j-5:j+4),2) nanmean(cmip6_ssp585_ts_vari2(:,j-5:j+4),2)];
% -- emergent constraint (following Brient and Schneider 2016: https://github.com/tapios/Emergent-constraint-regression)
[tmp_vari2_sort,tmp_vari2_sortidx] = sort(tmp_vari2);
tmp_sort = tmp(tmp_vari2_sortidx);
% -- OLS
[p, yhat, ci] = polypredci(tmp_vari2_sort,tmp_sort, 1, 0.95);

% -- York regression for regression of uncertain X and Y (and perhaps correlated errors) --
if strcmp(vari2,'siconc')==1
  r = 0*ones(1,length(tmp_vari2_sort)); % correlation coefficient between errors in X and Y
  std_predictor   = repmat(cmip6_piControl_ts_vari_std,[1,length(scen_cmip6)]);
  std_predictand  = std_predictor*0 + 0.56;
else
  r = repmat(cmip6_piControl_ts_vari_vari2_cor,[1,length(scen_cmip6)]); % internal variability correlation
  std_predictor   = repmat(cmip6_piControl_ts_vari_std,[1,length(scen_cmip6)]);
  std_predictand  = repmat(cmip6_piControl_ts_vari2_std,[1,length(scen_cmip6)]);
end
[a_york, b_york, sigma_a_york, sigma_b_york, b_save] = york_fit(tmp(:)', tmp_vari2(:)', std_predictor, std_predictand, r);
x_vals = sort(tmp);
yhat_york_vari2 = a_york + b_york * sort(tmp);
yhat_york_vari2_res = sort(tmp_vari2) - yhat_york_vari2;
yhat_york_vari2_sigma = nanstd(yhat_york_vari2_res);
% --
nboot       = 1e3;                % size of bootstrap sample;
nmodels     = length(models_cmip6);
tmp_constr      = [ssp119_loc2_ts_moore(j,2); ssp126_loc2_ts_moore(j,2); ssp245_loc2_ts_moore(j,2); ssp370_loc2_ts_moore(j,2); ssp585_loc2_ts_moore(j,2)];
std_tmp_constr  = sqrt((((ssp585_loc2_ts_moore(j,3)-ssp119_loc2_ts_moore(j,1))/3.29).^2)');
tmp_constr_mean = ssp585_loc2_ts_moore(j,3) - ((ssp585_loc2_ts_moore(j,3)-ssp119_loc2_ts_moore(j,1))/2);
tmp_constr    = normrnd(tmp_constr_mean, std_tmp_constr, nboot, 1);
% -- just for sanity check
std_tmp       = sqrt((((prctile(cmip6_ssp585_ts_em(:,j),95)-prctile(cmip6_ssp119_ts_em(:,j),5))/3.29).^2)');
tmp_mean      = prctile(cmip6_ssp585_ts_em(:,j),95) - ((prctile(cmip6_ssp585_ts_em(:,j),95)-prctile(cmip6_ssp119_ts_em(:,j),5))/2);
tmp2          = normrnd(tmp_mean, std_tmp, nboot, 1);
tmp2          = normrnd(nanmean(tmp(:)), nanstd(tmp(:)), nboot, 1);
'------------------------------'
[prctile(cmip6_ssp119_ts_em(:,j),5) prctile(cmip6_ssp585_ts_em(:,j),95)]
prctile(tmp2,[5 95])
'------------------------------'
[ssp119_loc2_ts_moore(j,1) ssp585_loc2_ts_moore(j,3)]
prctile(tmp_constr,[5 95])
'------------------------------'

x_vals_constr = sort(tmp_constr);
clear('yhat_vari2_nboot','yhat_vari2_nboot_constr',...
'yhat_york_vari2_nboot','yhat_york_vari2_nboot_constr')
yhat_vari2_nboot              = NaN(nboot,1);
yhat_vari2_nboot_constr       = NaN(nboot,1);
yhat_york_vari2_nboot         = NaN(nboot,length(tmp));
yhat_york_vari2_nboot_constr  = NaN(nboot,nboot);
for i = 1:nboot
  % i
  ibm     = randsample(nmodels, nmodels, 1); % random selection of models
  if length(scen_cmip6)==5
    ib      = [ibm; ibm+1*nmodels; ibm+2*nmodels; ibm+3*nmodels; ibm+4*nmodels]; % make sure the same models are selected across scenarios
  else
    ib      = [ibm; ibm+1*nmodels; ibm+2*nmodels; ibm+3*nmodels]; % make sure the same models are selected across scenarios
  end
  % -- OLS regression:
  [B, BINT, R, RINT, STATS]     = regress(tmp_vari2(ib), [ones(size(ib)) tmp(ib)]);
  sigma                         = sqrt(STATS(4));
  yhat_vari2_nboot(i)           = B(1)+B(2)*nanmean(tmp(ib)) + sigma*randn(1);  % prediction from mean model response
  yhat_vari2_nboot_constr(i)    = B(1)+B(2)*nanmean(tmp_constr(ib)) + sigma*randn(1);  % prediction from range of constrained model responses
  % -- York regression:
  if strcmp(vari2,'siconc')==0
    r = repmat(cmip6_piControl_ts_vari_vari2_cor(ibm),[1,length(scen_cmip6)]); % internal variability correlation
    std_predictor   = repmat(cmip6_piControl_ts_vari_std(ibm),[1,length(scen_cmip6)]);
    std_predictand  = repmat(cmip6_piControl_ts_vari2_std(ibm),[1,length(scen_cmip6)]);
  end
  [a_york, b_york, sigma_a_york, sigma_b_york, b_save] = york_fit(tmp(ib)', tmp_vari2(ib)', std_predictor, std_predictand, r);
  % -- predict
  yhat_york_vari2_nboot(i,:)        = a_york + b_york * tmp(:);
  yhat_york_vari2_nboot_constr(i,:) = a_york + b_york * tmp_constr;
end
% -- make sure sea ice is not less than 0:
if strcmp(vari2,'siconc')==1
  yhat_york_vari2_nboot(yhat_york_vari2_nboot<0) = 0.0001;
  yhat_york_vari2_nboot_constr(yhat_york_vari2_nboot_constr<0) = 0.0001;
end
% --
x = -20:50;
range = max(tmp_vari2)-min(tmp_vari2);
x = min(tmp_vari2)-0.1*range:max(tmp_vari2)+0.1*range;
% -- plot decadal mean summer tas
h0 = plot(tmp_vari2(:),tmp(:),'bo')%,'MarkerFaceColor','b');
h3 = plot(yhat_york_vari2, sort(tmp),'b') % York regression
[ci5,idx5]    = sort(prctile(yhat_york_vari2_nboot,5));
[ci95,idx95]  = sort(prctile(yhat_york_vari2_nboot,95));
plot(ci5, tmp(idx5),'b--') % York regression
plot(ci95, tmp(idx95),'b--') % York regression
plot([-0.5 1.1*max(tmp_vari2)],[0 0],'k')
% --
sf = 3;
if strcmp(region,'global')==1 || strcmp(region,'Central_Europe')==1
  sf = .5;
end
% -- PDF of unconstrained projection (2nd variable):
[f0,x0,u] = ksdensity(nanmean(yhat_york_vari2_nboot,1));
if strcmp(vari2,'siconc')==1
  idx = find(x0>=0);
  h1 = line(x0(idx),f0(idx)*sf,'Color',[.5 .5 1], 'LineWidth',2);
else
  h1 = line(x0,f0*sf,'Color',[.5 .5 1], 'LineWidth',2);
end
% -- PDF of constrained projection (2nd variable):
[f0,x0,u] = ksdensity(nanmean(yhat_york_vari2_nboot_constr,1));%,'Support','positive');
if strcmp(vari2,'siconc')==1
  idx = find(x0>=0);
  h2 = line(x0(idx),f0(idx)*sf,'Color',[1 .5 .5], 'LineWidth',2);
else
  h2 = line(x0,f0*sf,'Color',[1 .5 .5], 'LineWidth',2);
end
% -- horizontal bars
plot([prctile(nanmean(yhat_york_vari2_nboot,1),5) prctile(nanmean(yhat_york_vari2_nboot,1),95)],[ylim(1)+abs(ylim(1)*0.6) ylim(1)+abs(ylim(1)*0.6)],'Color',[.5 .5 1],'LineWidth',8)
plot([prctile(yhat_york_vari2_nboot_constr(:),5) prctile(yhat_york_vari2_nboot_constr(:),95)],[ylim(1)+abs(ylim(1)*0.2) ylim(1)+abs(ylim(1)*0.2)],'Color',[1 .5 .5],'LineWidth',8)
legend([h0 h3 h1 h2],['CMIP6 ' num2str(end_year-4) '-' num2str(end_year+5)],['r = ' num2str(round(corr(tmp,tmp_vari2),2))],'Unconstrained','Constrained','location','best')
% legend boxoff
xlim = [min(tmp_vari2)-0.1*(max(tmp_vari2)-min(tmp_vari2)) 1.2*max(tmp_vari2)];
set(gca,'YLim',ylim,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))],'Layer','top','XLim',xlim)
hline([ssp585_loc2_ts_moore(j,3), ssp119_loc2_ts_moore(j,1)],'k:')
vline([prctile(yhat_york_vari2_nboot_constr(:),5) prctile(yhat_york_vari2_nboot_constr(:),95)],'k:')
if refende <= 2020
  if strcmp(vari2,'tas')==1 && strcmp(seas2,'JJA')==1
    vline(0,'r') % 2001-2020
    vline(obs2(end),'r')
    text(0-.1,ylim(2)/3,['Average summer ' num2str(refstart) '-' num2str(refende)],'rotation',90,'color','r','FontSize',8)
    text(obs2(end)-.1,ylim(2)/2,'Summer 2022','rotation',90,'color','r','FontSize',8)
  end
end
box on
if strcmp(vari2,'pr')==1
  xlabel('Precipitation change (%)')
elseif strcmp(vari2,'tas')==1 && strcmp(seas2,'JJA')==1
  xlabel('JJA temperature change (\circC)')
elseif strcmp(vari2,'siconc')==1 && strcmp(seas2,'S')==1
  xlabel('September sea ice area (million km^2)')
end
% --
yyaxis right
set(gca,'YLim',ylim+offset,'YTick',[ceil(ylim(1)):yincr:round(ylim(end))*2],'XLim',xlim)
ylabel(['relative to preindustrial (\circC)'])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
% --
if refende <= 2020 && strcmp(vari2,'tas')==1
  handle = subplot(1,3,3);
  pos = get(handle,'position');
  b=axes('Position',[pos(1) 0. pos(3) 1e-12]);
  offset2 = abs(nanmean(obs2(1880-time_obs2(1)+1:1910-time_obs2(1))));
  set(b,'Color','none','Xlim',xlim+offset2);
  xlabel(b,'relative to preindustrial (\circC)')
end


tightfig

% return
set(gcf,'PaperPositionMode','auto');
set(gcf,'renderer','Painters')
% fileo = ['~/Dropbox/publication/lehner22_reduced_uncertainty/fig/hawkins_plots_cmip6_' region '_' vari '_' seas '_' num2str(wl) 'yr_ref' num2str(refstart) '-' num2str(refende) '_reduced_uncertainty_end' num2str(end_year) '_pie_scale'];
fileo = ['~/Dropbox/publication/lehner22_reduced_uncertainty/fig/lehner23aguadv_fig1'];
print('-r300','-loose', '-depsc', ['' fileo '.eps'])
save2pdf(['' fileo '.pdf'])
