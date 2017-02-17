function find_zi(data,recal)
clc; close all


if ~recal && exist('zi.mat','file')
    load('zi.mat')
end

%% download data from U wyoming servers
% download data
set(0,'defaulttextinterpreter','latex')
set(0,'defaultlinelinewidth',2)
stationNames = {'Belize','Cancun'};
close all
if 0
    stationIDs = {'78583', '76595'};
    saveFolder = 'C:\Users\Derek\Dropbox (Personal)\Paper 3\soundingData\';
    
    % Belize Airport: 78583
    % Cancun: 76595
    station='78583';
    region='naconf';
    
    % save data 2013
    years=2013;
    months=11:12;
    downloadRWS_URL(years,months,station,region,saveFolder)
    
    % save data 2014 - 2015
    years=2014:2015;
    months=1:5;
    downloadRWS_URL(years,months,station,region,saveFolder)
    
    % save data 2016
    years=2016;
    months=1:4;
    downloadRWS_URL(years,months,station,region,saveFolder)
    
    
    for ii = 1:2
        
        % find station
        station = stationIDs{ii};
        
        % find all files for station
        allFiles = dir([saveFolder,'*',station,'*']);
        
        % initialize output
        time = [];
        M = nan(100,11,[]);
        
        % load monthly files
        for jj = 1:numel(allFiles)
            [timeLocal, Mlocal] = read_RWS_html2([saveFolder,allFiles(jj).name],station);
            time(end+1:end+size(timeLocal,2),1) = timeLocal;
            M(:,:,end+1:end+size(Mlocal,3)) = Mlocal;
        end
        
        % store in structure
        soundingData.(stationNames{ii}).M = M;
        soundingData.(stationNames{ii}).time = time;
        
    end
    
    save('soundingData','soundingData')
end
%% Estimate zi
% load sounding data
load('soundingData')
Ribc = 0.25;

if ~exist('zi','var')
    
    % see http://weather.uwyo.edu/upperair/columns.html for column descriptions
    header = {'PRES'   'HGHT'   'TEMP'   'DWPT'   'RELH'   'MIXR'   'DRCT'   'SKNT'   'THTA'   'THTE'   'THTV'};
    
    % initialize variables
    zi.(stationNames{1}).z00.zi = [];
    zi.(stationNames{1}).z00.time = [];
    zi.(stationNames{1}).z12.zi = [];
    zi.(stationNames{1}).z12.time = [];
    zi.(stationNames{2}).z00.zi = [];
    zi.(stationNames{2}).z00.time = [];
    zi.(stationNames{2}).z12.zi = [];
    zi.(stationNames{2}).z12.time = [];
    
    numPlots = 1;
    for ii = 1:1
        
        for jj = 1:size(soundingData.(stationNames{ii}).time,1)
            try
                
                % find local M and local time
                localM = soundingData.(stationNames{ii}).M(:,:,jj);
                localTime = soundingData.(stationNames{ii}).time(jj);
                
                % only looking at April to August
                [~,mm,~] = datevec(floor(localTime));
                if mm < 4 || mm > 8
                    continue
                end
                
                if localTime - floor(localTime) > 0.1 % only care about 0z
                    continue
                end
                
                % find rows up to 3 km
                rows = find(localM(:,2) < 3000);
                
                if isempty(rows)
                    continue
                end
                
                % find z
                z = localM(rows,2);
                
                % find T
                T = localM(rows,3)+273.15;
                
                % find theta_v
                theta_v = localM(rows,11);
                
                % check for unstable layers
                if theta_v(1) < theta_v(2)
                    continue
                end
                
                % find theta_v0
                theta_v0 = theta_v(1);
                
                % interpolate theta_v
                theta_v = interp1(z,theta_v,z(1):z(end))';
                
                % smooth theta_v
                theta_v = smooth(theta_v,300,'loess');
                
                % find speed
                spd = localM(rows,8) * 0.5144; % knots to m/s
                
                % interpolate speed
                spd = interp1(z,spd,z(1):z(end))';
                
                % smooth speed
                spd = smooth(spd,300,'loess');
                
                % find RH
                RH = localM(rows,5);
                
                % interpolate RH
                RH = interp1(z,RH,z(1):z(end))';
                
                % smooth RH
                RH = smooth(RH,300,'loess');
                
                % redefine z
                z = (z(1):z(end))';
                
                % find Rib
                Rib = (9.81/theta_v0)*(theta_v - theta_v0).*z ./ spd.^2;
                
                % find gradient
                dTdz = gradient(theta_v);
                
                % find zi
                ziRib = z(find(Rib>Ribc,1,'first'));
                
                % find zi from gradient
                [~, index] = max(dTdz);
                ziMaxGrad = z(index);
                
                if ziRib < 400 || ziRib > 2000
                    continue
                end
                
                % plot data
                if 1 % mod(jj,30) < 2
                    fullFigure;
                    %                 subplot1(1,5,'Gap',[0 0])
                    subplot1(1,2,'Gap',[0 0])
                    changeFig(20,2,10)
                    %                 subplot1(1)
                    %                 h = plot(Rib,z,[Ribc Ribc],[z(1) z(end)],'k--');
                    %                 set(h(2),'linewidth',1)
                    %                 ylabel('$z$ (m)','interpreter','latex')
                    %                 xlabel('$Rib$')
                    %                 set(gca,'xtick',[min(xlim)+0.2*range(xlim) mean(xlim) min(xlim)+0.8*range(xlim)])
                    
                    %                 subplot1(2)
                    %                 plot(T,localM(rows,2),'bo-',localM(rows,11),localM(rows,2),'go-',theta_v,z,nanmean(theta_v),ziRib,'k+',nanmean(theta_v),ziMaxGrad,'bs')
                    %                 h = legend('$T$ raw','$\theta_v$ raw','$\theta_v$ smooth','$zi_{\mathrm{Rib}}$','$zi_{\mathrm{Max~Grad}}$','location','northwest');
                    %                 xlabel('K')
                    %                 set(h,'interpreter','latex')
                    %                 set(gca,'xtick',[min(xlim)+0.2*range(xlim) mean(xlim) min(xlim)+0.8*range(xlim)])
                    %
                    %                 subplot1(3)
                    %                 plot(dTdz,z)
                    %                 title([stationNames{ii},' - ',datestr(localTime,'dd mmm yy - hhz')])
                    %                 xlabel('$\frac{\partial \theta_v}{\partial z}$','interpreter','latex')
                    %                 set(gca,'xtick',[min(xlim)+0.2*range(xlim) mean(xlim) min(xlim)+0.8*range(xlim)])
                    %
                    %                 subplot1(4)
                    %                 plot(localM(rows,5),localM(rows,2),'ko',RH,z)
                    %                 xlabel('$RH$ (\%)','interpreter','latex')
                    %                 set(gca,'xtick',[min(xlim)+0.2*range(xlim) mean(xlim) min(xlim)+0.8*range(xlim)])
                    %
                    %                 subplot1(5)
                    %                 plot(localM(rows,8)*0.5144,localM(rows,2),'ko',spd,z)
                    %                 xlabel('$\frac{\partial \theta_v}{\partial z}$','interpreter','latex')
                    %                 set(gca,'xtick',[min(xlim)+0.2*range(xlim) mean(xlim) min(xlim)+0.8*range(xlim)])
                    subplot1(1)
                    plot(T,localM(rows,2),'bo-',nanmean(T),ziRib,'k+',nanmean(T),ziMaxGrad,'bs')
                    h = legend('$T$ raw','$zi_{\mathrm{Rib}}$','$zi_{\mathrm{Max~Grad}}$','location','north');
                    title([stationNames{ii},' - ',datestr(localTime,'dd mmm yy - hhz')])
                    xlabel('K')
                    set(h,'interpreter','latex')
                    set(gca,'xtick',[min(xlim)+0.2*range(xlim) mean(xlim) min(xlim)+0.8*range(xlim)])
                    
                    subplot1(2)
                    plot(localM(rows,11),localM(rows,2),'go-',theta_v,z,nanmean(theta_v),ziRib,'k+',nanmean(theta_v),ziMaxGrad,'bs')
                    h = legend('$\theta_v$ raw','$\theta_v$ smooth','$zi_{\mathrm{Rib}}$','$zi_{\mathrm{Max~Grad}}$','location','northwest');
                    xlabel('K')
                    set(h,'interpreter','latex')
                    set(gca,'xtick',[min(xlim)+0.2*range(xlim) mean(xlim) min(xlim)+0.8*range(xlim)])
                    
                    % check for correction
                    correctZiFlag = input('Correct ziRib? 0 for no, 1 for yes:  ');
                    if correctZiFlag
                        display('Locat zi on right panel of figure')
                        ziObs = ginput(1);
                        ziObs = ziObs(2);
                    else
                        ziObs = [];
                    end
                    
                    if ~isempty(ziObs)
                        cla
                        plot(localM(rows,11),localM(rows,2),'go-',theta_v,z,nanmean(theta_v),ziRib,'k+',nanmean(theta_v),ziMaxGrad,'bs',nanmean(theta_v),ziObs,'rx')
                        h = legend('$\theta_v$ raw','$\theta_v$ smooth','$zi_{\mathrm{Rib}}$','$zi_{\mathrm{Max~Grad}}$','$zi_{\mathrm{Manual}}$','location','northwest');
                        xlabel('K')
                        set(h,'interpreter','latex')
                        set(gca,'xtick',[min(xlim)+0.2*range(xlim) mean(xlim) min(xlim)+0.8*range(xlim)])
                        display('zi manual will be used')
                        ziRib = ziObs;
                    end
                    
                    
                    % export fig
                    set(gcf,'color','white')
                    export_fig([stationNames{ii},'_',datestr(localTime,'hhz')],'-pdf','-append')
                end
                
                if localTime - floor(localTime) > 0.1 % 4:00 AM
                    zi.(stationNames{ii}).z12.zi(end+1) = ziRib;
                    zi.(stationNames{ii}).z12.time(end+1) = localTime;
                else % 6:00 PM
                    zi.(stationNames{ii}).z00.zi(end+1) = ziRib;
                    zi.(stationNames{ii}).z00.time(end+1) = localTime;
                end
                
                %numPlots = numPlots+1;
                %if numPlots > 6
                %    break
                %end
                
            catch err
                err.message
            end
        end
    end
    
    saveFlag = input('Save data true or false:  ');
    if saveFlag
        % save zi data
        save('zi','zi')
    end
    
    fullFigure;
    plot(zi.Belize.z00.time,zi.Belize.z00.zi,'rsq')
    hold on
    plot(zi.Belize.z12.time,zi.Belize.z12.zi,'bsq')
    plot(zi.Cancun.z00.time,zi.Cancun.z00.zi,'r*')
    plot(zi.Cancun.z12.time,zi.Cancun.z12.zi,'b*')
    changeFig(24,3)
    legend('Belize 0z','Belize 12z', 'Cancun 0z', 'Cancun 12z')
    axis([datenum(2013,12,1) datenum(2016,5,1) 400 2000])
    datetick('x','mmm yy')
    set(gcf,'color','white')
    export_fig('ziEstimates','-pdf')
end
%% find inputs for tke model

% find possible dates
dates = zi.Belize.z00.time;

% iterate through sites
for ii = 1:2
    
    if ii == 1; site = 'beach'; else site = 'substation'; end
    
    % iterate through dates and create structures
    for jj = 1:numel(dates)
               
        % find local date
        dateLocal = zi.Belize.z00.time(jj);
        
        onshoreFlag = false;
        offshoreFlag = false;
        
        % check for onshore flow
        if data.(site).dailyWindMode(data.allDays == dateLocal) > data.onshoreDir - data.dtheta && data.(site).dailyWindMode(data.allDays == dateLocal) < data.onshoreDir + data.dtheta
            display(sprintf('%s: %s onshoreFLow:  %g',datestr(dateLocal,'dd mmm yy'),site,data.(site).dailyWindMode(data.allDays == dateLocal)))
            onshoreFlag = true;
        elseif data.(site).dailyWindMode(data.allDays == dateLocal) > data.offshoreDir - data.dtheta && data.(site).dailyWindMode(data.allDays == dateLocal) < data.offshoreDir + data.dtheta
            offshoreFlag = true;
            continue
            display(sprintf('%s: %s offshoreFLow:  %g',datestr(dateLocal,'dd mmm yy'),site,data.(site).dailyWindMode(data.allDays == dateLocal)))
        else
            %display(sprintf('Neither onshore nor offshore:  %g',data.(site).dailyWindMode(data.allDays == dateLocal)))
            continue
        end
            
        % zi
        ziLocal = zi.Belize.z00.zi(jj);
        
        % find day rows
        rows = find(dateLocal == floor(data.tLST)) + 1;
        
        % wind
        windSpdCol10 = strcmp(data.(site).min10Header,'cup10m_Avg');
        windSpdCol20 = strcmp(data.(site).min10Header,'cup20m_Avg');
        windSpdCol40 = strcmp(data.(site).min10Header,'cup40m_Avg');
        WS10 = data.(site).min10(rows,windSpdCol10);
        WS20 = data.(site).min10(rows,windSpdCol10);
        WS40 = data.(site).min10(rows,windSpdCol10);
        
        % sensible heat flux
        kinH = data.(site).H(rows,5);
        
        % wind
        SWcol = strcmp(data.(site).min10Header,'solar_Avg');
        SW = data.(site).min10(rows,SWcol);
        
        % temperature
        TCol = strcmp(data.(site).min10Header,'T10m_Avg');
        T10m = data.(site).min10(rows,TCol) + 273.15;
        
        % check for nans
         %sum(isnan(WS))/numel(dateLocal)>0.1 || sum(isnan(kinH))/numel(dateLocal)>0.1 || sum(isnan(T10m))/numel(dateLocal)>0.1;
        if sum(isnan(WS10))/numel(dateLocal)>0.1 || sum(isnan(kinH))/numel(dateLocal)>0.1 || sum(isnan(T10m))/numel(dateLocal)>0.1
            continue
        end
        
        % zi
        tkeModelInputs.(site).(datestr(dates(jj),'mmm_dd_yy')).zi = ziLocal;
        
        % time
        tkeModelInputs.(site).(datestr(dates(jj),'mmm_dd_yy')).tLST = data.tLST(rows);
        
        % wind
        tkeModelInputs.(site).(datestr(dates(jj),'mmm_dd_yy')).WS10 = WS10;
        tkeModelInputs.(site).(datestr(dates(jj),'mmm_dd_yy')).WS20 = WS20;
        tkeModelInputs.(site).(datestr(dates(jj),'mmm_dd_yy')).WS40 = WS40;
        
        % sensible heat flux
        tkeModelInputs.(site).(datestr(dates(jj),'mmm_dd_yy')).kinH = kinH;
        
        % temperature
        tkeModelInputs.(site).(datestr(dates(jj),'mmm_dd_yy')).T10m = T10m;
        
        % daily wind mode
        tkeModelInputs.(site).(datestr(dates(jj),'mmm_dd_yy')).windMode = data.(site).dailyWindMode(data.allDays == dateLocal);
        
        % wind flags
        tkeModelInputs.(site).(datestr(dates(jj),'mmm_dd_yy')).onshoreFlag = onshoreFlag;
        tkeModelInputs.(site).(datestr(dates(jj),'mmm_dd_yy')).offshoreFlag = offshoreFlag;
        
        % solar 
        tkeModelInputs.(site).(datestr(dates(jj),'mmm_dd_yy')).solar = SW;
        tkeModelInputs.(site).(datestr(dates(jj),'mmm_dd_yy')).solar = SW;
    end
end

save('tkeModelInputs','tkeModelInputs')

end



