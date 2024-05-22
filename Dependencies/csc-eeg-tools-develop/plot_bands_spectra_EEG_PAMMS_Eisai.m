function [psd_average]=plot_bands_spectra_EEG_PAMMS_Eisai(EEG)


EEG.psd.Bands.SWA     = find(EEG.psd.Hzbins > 0.35   & EEG.psd.Hzbins    <= 4);
EEG.psd.Bands.Theta   = find(EEG.psd.Hzbins > 4   & EEG.psd.Hzbins    <= 8);
EEG.psd.Bands.Alpha   = find(EEG.psd.Hzbins > 8   & EEG.psd.Hzbins    <= 12);
EEG.psd.Bands.Sigma   = find(EEG.psd.Hzbins > 12  & EEG.psd.Hzbins    <= 15);
EEG.psd.Bands.Beta    = find(EEG.psd.Hzbins > 15  & EEG.psd.Hzbins    <= 25);
EEG.psd.Bands.Gamma   = find(EEG.psd.Hzbins > 25  & EEG.psd.Hzbins    <= 35);
%EEG.psd.Bands.All     = find(EEG.psd.Hzbins > 1   & EEG.psd.Hzbins    <= 40);
Bands               = fieldnames(EEG.psd.Bands);
numbands    = length(Bands);
%numchannels = size(EEG.psd.data,1);
numtrials   = size(EEG.psd.data,4);


psd_average     = squeeze(nanmedian(EEG.psd.data,3)); % average across epochs

for n = 1:numtrials
topofig = figure('position',[20 50 500 150*numbands]);
subplot('position',[.1 .96 .8 .02])
figlabel = strrep(EEG.setname,'_',' ');
text(.5,.5,figlabel,'FontSize',12,'fontweight','bold','horizontalalignment','center');
%text(.1,.5,['Trial ',num2str(n)],'FontSize',12,'fontweight','bold','horizontalalignment','center');
colormap(jet);
axis off
outlier_tot = [];

    for b=1:numbands
        figure(topofig);
        
        subplot('position',[.15 (1-.07)*(numbands-b)/numbands .7 (1-.07)/numbands]);
        data = squeeze(nanmean(psd_average(:,EEG.psd.Bands.(Bands{b}),n),2));
        % this determines outliers - could be put in as an input
        outlier = find(bsxfun(@gt, abs(bsxfun(@minus, data,mean(data))), 3*std(data)));
        
        topoplot(data,EEG.chanlocs,'shading','interp','style','map','maplimits','maxmin','whitebk','on','electrodes','on');
        electrodes.x=get(findobj(gca,'Marker','.'),'XData'); electrodes.y=get(findobj(gca,'Marker','.'),'YData'); electrodes.z=get(findobj(gca,'Marker','.'),'ZData');
        delete(findobj(gca,'Marker','.'));
        for oi = 1:length(outlier)
            outch = outlier(oi);
            text(electrodes.x(outch),electrodes.y(outch),EEG.chanlocs(outch).labels,'Fontsize',8);
        end
        set(gca,'Xlim',[-.55 .55]); set(gca,'Ylim',[-.59 .59]);
        text(-.8,0,(Bands{b}),'FontSize',18,'rotation',90,'HorizontalAlignment','center');
        h=colorbar;
        set(h,'position',[.68 (1-.05)*(numbands-b)/numbands .05 (1-.2)/numbands]);
        outlier_tot = unique([outlier_tot;outlier]);  
    end;
    
    set(topofig,'color','w','paperpositionmode','auto');
    print(topofig,'-dpng','-r500',[EEG.filepath,filesep,EEG.setname,'_topo.png']);
   
  
    globalspectralfig = figure('position',[550 200 550 900]);
    h=semilogy(EEG.psd.Hzbins,psd_average);
    title(EEG.setname,'FontSize',24,'fontweight','bold','horizontalalignment','center');
    colormap(jet);
    
    for c = 1:EEG.nbchan
        hold on;
        text(EEG.psd.Hzbins(end),psd_average(c,end),EEG.chanlocs(c).labels);
        xlim([0 EEG.psd.Hzbins(end)+5]);
    end
    
    for oi = 1:length(outlier_tot)
            outch = outlier_tot(oi);
            set(h(outch),'linewidth',3)
    end
    
    set(globalspectralfig,'color','w','paperpositionmode','auto')
    
    
    print(globalspectralfig,'-dpng','-r500',[EEG.filepath,filesep,EEG.setname,'_spectra.png']);
    savefig(globalspectralfig,[EEG.filepath,filesep,EEG.setname,'_spectra']);
   
end
