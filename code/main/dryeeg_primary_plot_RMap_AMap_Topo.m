%% DRy EEG Figures
% Author: Bahne Bahners
clear variables
% Repository-relative setup (original source remains unchanged).
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repoRoot, 'functions'), '-begin');
dryeeg_setup(repoRoot);
% Machine-specific addpath removed; see config/local_paths.m
% Machine-specific addpath removed; see config/local_paths.m
% Machine-specific addpath removed; see config/local_paths.m
cohmappath=[dryeeg_data_path(repoRoot, 'EPmaps_Discovery'), filesep]; % connectivity map path
rmappath=dryeeg_result_path(repoRoot, 'figures_RMap_AMap_Topo'); % results path
% Make sure resultpath exists
if ~exist([rmappath])
    mkdir([rmappath]);
end
% Machine-specific addpath removed; see config/local_paths.m
load(dryeeg_data_path(repoRoot, 'mandrillcolormap.mat')); 
load(dryeeg_data_path(repoRoot, 'timevec.mat')); 
load(dryeeg_data_path(repoRoot, 'channels.mat')); 
corrtype='Spearman'; % use which type of correlation metric - for fMRI could also do Pearson

%% Load Regressor Table
tab=readtable(dryeeg_data_path(repoRoot, 'metadata/UPDRS_dryEEG_ImprovementsOmitNanUp.xlsx'));
tab([27],:)=[]; % P027 (left out for external validation)
allpts=1:size(tab,1); % all patients vector (cave: only used for first loop), then hemispheres (allconds)

%% Analysis time window:
ini=0.01; % beginning of time window in seconds
fin=0.2 ; % end of time window in seconds
win=[min(find(time>=ini)):min(find(time>=fin))];
Tmap=repmat(time(win),[32,1]); %time map
time=time(win); % time based on window defined

%% Prepare Channel Order & Channel Flip
channels(:,33:35)=[];
for ii=1:size(channels,2)
    channlabels{ii}=channels(ii).Name; % create cell array with EEG channel labels
end
chanord= bb_chanord(channlabels);% now get a new channelorder for the plots (Frontal to Occipital)
chanlist_new=char(channlabels{chanord}); % now reorder the channellist
% new channellabels are assigned to essentially flip the EEG data
chanside={'right','right','right','right','right','right','right','right','right','right','left','left','z','right','z','right','left','left','z','right','right','z','left','left','left','left','left','left','left','left','left','left'};
chanflip={'P7',	'T7', 'CP5',	'FC5',	'F7',	'F3',	'C3',	'P3',	'AF3',	'Fp1',	'Fp2',	'AF4',	'Fz',	'FC1',	'Cz',	'CP1',	'PO4',	'O2',   'Oz',	'O1',	'PO3',	'Pz',	'CP2',	'FC2',	'P4',	'C4',	'F4',	'F8',	'FC6',	'CP6',	'T8',	'P8'};
[~, flip] = ismember(chanflip,channlabels); % Find the corresponding indices for channel reordering
ridx=contains(chanside,'right');
lidx=contains(chanside,'left');

%% 0. create a matrix of EEG Evoked Potential Files
for pt=allpts % iterate through patients - this loop will just give us a cell with entries that map to each patients structural connectivity map (seeding from the VTAs).
    list=dir([cohmappath]);
    files=char(list.name);
    clear fn
    fn=files(contains(string(files),tab.SubID{pt}),:);
    left=load([list(pt).folder,'/',regexprep(fn(contains(string(fn),'le'),:),' ','')]);
    right=load([list(pt).folder,'/',regexprep(fn(contains(string(fn),'ri'),:),' ','')]);
    L=left.m(:,win);
    %L(ridx,:)=NaN;
    Ri=right.m(:,win);
    %Ri(lidx,:)=NaN;
    patConnectivityFiles{pt,1}=Ri(flip,:);
    patConnectivityFiles{size(allpts,2)+pt,1}=L;
    Regressor{pt,1}=tab.DeltaPercL(pt); % percentage improvement right hemiscores
    Regressor{size(allpts,2)+pt,1}=tab.DeltaPercR(pt); % percentage improvement left hemiscores
    Groupcell{pt}='1';

end

Regressor=cell2mat(Regressor);
channelvec=1:size(patConnectivityFiles{1},1);
allconds=1:size(Regressor,1);
%RegressorHat=NaN(size(allconds,2),1);
    %% 1. generate R-Maps
    % usage: Rmap=bb_Rmap_ep(fingerprints,regressor,output,corrtype)
    Rmap=bb_Rmap_ep(patConnectivityFiles(allconds),Regressor(allconds),...
        [rmappath,'R_',sprintf('%02.0f',0),'.mat'],corrtype);

   h= ea_corrplot(Regressor,Regressor)
    

Chmap=repmat(channelvec,[size(patConnectivityFiles{1},2),1])';

% Time windows of interest
% Time window 1
init=0.02;
fint=0.03;
win1=[min(find(time>=init)):min(find(time>=fint))];
% Time window 2
init=0.05;
fint=0.06;
win2=[min(find(time>=init)):min(find(time>=fint))];
% Time window 3
init=0.075;
fint=0.085;
win3=[min(find(time>=init)):min(find(time>=fint))];
% Time window 4
init=0.09;
fint=0.12;
win4=[min(find(time>=init)):min(find(time>=fint))];

timwin={win1,win2,win3,win4};
timwinlab={'   20 to 30 ms', '   50 to 60 ms','   75 to 85 ms','   90 to 120 ms'};
intp=1000;
% Topoplots
halfidx=contains(chanside,'left')|contains(chanside,'z') ;
h2=figure;
for iwin=1:size(timwin,2)
    a= subplot(2,4,iwin)

    [h,ch_x,ch_y]=plot_topography(channlabels,(mean(Rmap(1:32,timwin{iwin})')),false(1),'10-20',false(1),false(1),intp);%,'false')%,1,1,0)
    hold on
    scatter(ch_x(halfidx,:), ch_y(halfidx,:), 10,'white', 'LineWidth',1);
    
    fg=gcf;
    fg.Colormap=cmap;
    title(timwinlab{iwin},'FontSize',10);
    cb = colorbar();
    clim([-0.4 0.4]);
    cb.Visible='off';
    ax=gca;
    ax.TickDir="out";
    ax.Children(1).LineWidth=1;
    ax.Children(2).LineWidth=1;
    ax.Children(3).LineWidth=1;
    ax.Children(4).LineWidth=1;
    a.Position (2)=a.Position(2)+0.05%1
    %     if iwin==1
    %     a.Position (1)=a.Position(1)-0.06
    %     elseif iwin==2
    % a.Position (1)=a.Position(1)-0.13
    %     elseif iwin==3
    a.Position (1)=a.Position(1)-0.06
    %end
end
a=subplot(2,4,5:8) % R-Matrix

Rmap(~halfidx,:)=NaN;
Rmap=Rmap(chanord,:);
chanlist_new(sum(isnan(Rmap),2)>1,:)=[];
channelvechalf=1:length(chanlist_new);
Rmap(sum(isnan(Rmap),2)>1,:)=[];


Chmap=repmat(channelvechalf,[size(patConnectivityFiles{1},2),1])';

imagesc(time,channelvechalf,Rmap)
colormap(cmap)
clim([min(min(Rmap)),max(max(Rmap))]);
xlabel('Time (s)','FontSize',16)
yticklabels(chanlist_new)
yticks(channelvechalf)
cb = colorbar();
ylabel(cb,'R','FontSize',10,'Rotation',0)
%cb.Label.Position(1) = 3;
cb.FontSize=14;
cb.FontWeight='bold';
h2.Color='white';
ax=gca;
ax.TickDir="out";
ax.XAxis.FontSize=14;
ax.XAxis.FontWeight="bold";
ax.FontWeight='bold';


%Position move %% position =[x_position y_position widht length] all are in some unit
h2.Position=h.Position;
%a.Position(4)=a.Position(4)+0.27;
a.Position(4)=a.Position(4)+0.18;
a.Position(2)=a.Position(2)+0.03 ;
%h2.Position(4)=h2.Position(4)+100
%exportgraphics(h2,[rmappath,'RmapTopo.png'],'Resolution',600)
saveas(h2,[rmappath,'RmapTopo.fig']);
saveas(h2,[rmappath,'RmapTopo.png']);

set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gca,'Box','off','Color','none')

export_fig([rmappath,'RmapTopo_transparent.png'], '-png', '-transparent');
% set(gcf, 'Color', 'none');
% set(gca, 'Color', 'none');
% export_fig([rmappath,'RmapTopo_transparent.tif'], '-tif', '-transparent','-r800');

chanlist_new=char(channlabels{chanord});

%% 5. Grand Average Time Series Plot 
% Redefine EP map time and reload maps
% Analysis time window:
load(dryeeg_data_path(repoRoot, 'timevec.mat')); 
ini=-0.03; % beginning of time window in seconds
fin=0.2 ;% end of time window in seconds
win=[min(find(time>=ini)):min(find(time>=fin))];
Tmap=repmat(time(win),[32,1]); %time map

time=time(win); % time based on window defined

% 0. create a new matrix of EEG Evoked Potential Files
for pt=allpts % iterate through patients - this loop will just give us a cell with entries that map to each patients structural connectivity map (seeding from the VTAs).
    list=dir([cohmappath]);
    files=char(list.name);
    clear fn
    fn=files(contains(string(files),tab.SubID{pt}),:);
    left=load([list(pt).folder,'/',regexprep(fn(contains(string(fn),'le'),:),' ','')]);
    right=load([list(pt).folder,'/',regexprep(fn(contains(string(fn),'ri'),:),' ','')]);
    L=left.m(:,win);
    R=right.m(flip,win);
    patConnectivityFiles{pt,1}=R;
    patConnectivityFiles{size(allpts,2)+pt,1}=L;
end
Amap=mean(cat(3,patConnectivityFiles{:}),3);
Smap=std(cat(3,patConnectivityFiles{:}),0,3)/ sqrt( size(patConnectivityFiles,1) );

% Load the corrected baseline-permutation mask produced by
% dryeeg_evoked_vs_baseline_permutation.m. Exact time-point matching keeps
% the overlay aligned with the plotted grand-average map.
permutationMask=false(size(Amap));
statsFile=fullfile(repoRoot,'results','baseline_statistics', ...
    'evoked_vs_baseline_statistics.mat');
if exist(statsFile,'file')
    statsLoaded=load(statsFile,'results');
    statsTime=round(statsLoaded.results.time(:)',10);
    figureTimeRounded=round(time(:)',10);
    [timeFound,statsLocation]=ismember(figureTimeRounded,statsTime);
    permutationMask(:,timeFound)= ...
        statsLoaded.results.significantMask(:,statsLocation(timeFound));
else
    warning(['Baseline-permutation statistics were not found. Run ' ...
        'dryeeg_evoked_vs_baseline_permutation.m before this script.']);
end
% Time windows of interest
% Time window 1
init=0.020;
fint=0.030;
win1=[min(find(time>=init)):min(find(time>=fint))];
% Time window 2
init=0.05;
fint=0.06;
win2=[min(find(time>=init)):min(find(time>=fint))];
% Time window 3
init=0.075;
fint=0.085;
win3=[min(find(time>=init)):min(find(time>=fint))];
% Time window 4
init=0.090;
fint=0.12;
win4=[min(find(time>=init)):min(find(time>=fint))];

timwin={win1,win2,win3,win4};
timwinlab={'   20 to 30 ms', '   50 to 60 ms','   75 to 85 ms','   90 to 120 ms'};
%% Plot starts here

% Define Color maps 
mand=mandrill(256);
mand=mand(end:-1:1,:);
channlabels_new=cellstr(regexprep(string(chanlist_new),' ',''));
%vcmap=viridis(size(Amap,1));
vcmap=mandrill(size(Amap,1)+10);
vcmap=vcmap(end:-1:1,:);
vcmap(19:23,:)=[];
%chanord=chanord(end:-1:1,:);

h4=figure;


for ch=1:size(Amap,1)
    %subplot(8,4,ch)
    plot(time,Amap(chanord(ch),:)*10e6,"Color",vcmap(ch,:))
    hold on
    conf= ciplot(Amap(chanord(ch),:)*10e6-Smap(chanord(ch),:)*10e6,Amap(chanord(ch),:)*10e6+Smap(chanord(ch),:)*10e6,time,vcmap(ch+5,:),0.2);
    conf.EdgeColor='none';
    %title(channlabels{chanord(ch)})
    ylim([-80,80]);
    xlim([ini,fin]);
end

xlabel('Time (s)','FontSize',14)
ylabel('Amplitude (ÂµV)','FontSize',14,'Rotation',90)
h4.Color='white';
ax=gca;
ax.TickDir="out";
ax.Box="off";
ax.XGrid="on";
ax.YGrid="on";
ax.FontSize=14;
ax.FontWeight="bold";
text(0.024, 67, 'P25','FontSize',14,'FontWeight','bold')
text(0.055, -40, 'N55','FontSize',14,'FontWeight','bold')
text(0.08, 15, 'P80','FontSize',14,'FontWeight','bold')
text(0.1, -20, 'N100','FontSize',14,'FontWeight','bold')

%set(gca, 'YDir','reverse')

ax.Colormap=vcmap;
c=colorbar;

c.Ticks=[0:1/31:1];
c.Direction="reverse";
c.TickLabels=channlabels_new;
c.TickDirection="out";
h4.Position=h.Position;
h4.Position(4)=h4.Position(4)+100;

%exportgraphics(h4,[rmappath,'GrandAverage.png'],'Resolution',600)
saveas(h4,[rmappath,'GrandAverage.fig']);

saveas(h4,[rmappath,'GrandAverage.png']);
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gca,'Box','off','Color','none')

export_fig([rmappath,'GrandAverage_transparent.png'], '-png', '-transparent');
%export_fig([rmappath,'GrandAverage_transparent.tif'], '-tif', '-transparent','-r800');


%% 4. Create Average Channels x Time Map
% Topoplots
h3=figure;
for iwin=1:size(timwin,2)
    a= subplot(2,4,iwin)
    plot_topography(channlabels,(mean(Amap(1:32,timwin{iwin})')*10e6),false(1),'10-20',false(1),false(1),intp)%,'false')%,1,1,0)
    fg=gcf;
    fg.Colormap=cmap;
    title(timwinlab{iwin},'FontSize',10);
    cb = colorbar();
    clim([-30 30])
    cb.Visible='off';
    ax=gca;
    ax.TickDir="out";
    ax.Children(1).LineWidth=1;
    ax.Children(2).LineWidth=1;
    ax.Children(3).LineWidth=1;
    ax.Children(4).LineWidth=1;
    a.Position(2)=a.Position(2)+0.1
    %     if iwin==1
    %     a.Position (1)=a.Position(1)-0.06
    %     elseif iwin==2
    % a.Position (1)=a.Position(1)-0.13
    %     elseif iwin==3
    a.Position (1)=a.Position(1)-0.06
    %end
end
a=subplot(2,4,5:8) % R-Matrix
imagesc(time(min(find(time>=ini)):min(find(time>=fin))),[1:32],Amap(chanord,:)*10e6)
colormap(cmap)
clim([-30,30]);
hold on
if any(permutationMask(:))
    contour(a,time,[1:32],double(permutationMask(chanord,:)), ...
        [0.5 0.5],'k','LineWidth',1.3);
end
xlabel('Time (s)','FontSize',10)
yticklabels(chanlist_new)
yticks([1:32])
cb = colorbar();
ylabel(cb,'Amplitude (ÂµV)','FontSize',10,'Rotation',90)
cb.Label.Position(1) = 3;
cb.FontSize=14;
cb.FontWeight='bold';
h3.Color='white';
ax=gca;
ax.TickDir="out";
ax.XAxis.FontSize=14;
ax.XAxis.FontWeight="bold";
ax.FontWeight='bold';

%Position move %% position =[x_position y_position widht length] all are in some unit
h3.Position=h.Position;
a.Position(4)=a.Position(4)+0.27 ;
a.Position(2)=a.Position(2)+0.015 ;
%sgtitle(regexprep(b(1:end-9),'_',' '),'FontWeight','bold')
h3.Position(4)=h3.Position(4)+100
%exportgraphics(h3,[rmappath,'AmapTopo.png'],'Resolution',600)
saveas(h3,[rmappath,'AmapTopo.fig']);
saveas(h3,[rmappath,'AmapTopo.png']);

set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gca,'Box','off','Color','none')

export_fig([rmappath,'AmapTopo_transparent.png'], '-png', '-transparent');

%export_fig([rmappath,'AmapTopo_transparent.tif'], '-tif', '-transparent','-r800');


