function [] = plotgraphs2()
% UNTITLED2 Summary of this function goes here
% Detailed explanation goes here

minvtopic = 2;
filename1 = ['finalbigrun6.mat'];
filename2 = ['finalbigrun6_std.mat'];
savefilepath = '/v/filer4b/v16q001/data/Ayans/Results_finalbigrun6/';

A = load(filename1);
B = load(filename2);
k2mat = [20:10:200];

for pind =1:10
        
    C = squeeze(A.savedstat(pind,:,:));
    C(:,3) = C(:,3);
    D = squeeze(B.savedstatstd(pind,:,:));
    %%C(:,2)
    %%D(:,2)

    p = [0.1:0.1:1]';
    
    rbaseline = 0.23*ones(size(k2mat));
    
    figure (1);
    
    %%errorbar(k2mat,C(:,1),D(:,1),'r.-', 'LineWidth',2, 'MarkerEdgeColor','r','MarkerSize',15); hold on;

    errorbar(k2mat,C(:,2),D(:,2),'g.-','LineWidth',2, 'MarkerEdgeColor','g','MarkerSize',15); hold on;
    errorbar(k2mat,C(:,3),D(:,3),'b.-','LineWidth',2, 'MarkerEdgeColor','b','MarkerSize',15); hold on;
    errorbar(k2mat,C(:,4),D(:,4),'k.-','LineWidth',2, 'MarkerEdgeColor','k','MarkerSize',15); hold on;
    errorbar(k2mat,C(:,5),D(:,5),'r.-','LineWidth',2, 'MarkerEdgeColor','r','MarkerSize',15); hold on;
    errorbar(k2mat,C(:,6),D(:,6),'m.-','LineWidth',2, 'MarkerEdgeColor','m','MarkerSize',15); hold on;
    plot(k2mat,rbaseline,'c.-', 'LineWidth',1.5); grid on, hold on;    
    
    t = xlabel('number of latent topics'); hold on;
    set(t, 'FontSize', 14);
    t = ylabel('over-all classification accuracy'); hold on;
    set(t, 'FontSize', 14);
    titlestr = ['Learning curve for p = ' num2str(p(pind))];
    t = title(titlestr);
    set(t, 'FontSize', 14);
    set(findobj('type','axes'),'fontsize',12);
    set(gca,'XLim',[0 220]);
    set(gca,'YLim',[0.1 1]);
    hold on;
    
    %%legend({'medLDA-OVA', 'DSLDA-NSLT', 'medLDA-MTL', 'DSLDA', 'DSLDA-OSST', 'random'}, 'Location','NorthWest'); 
    legend({'medLDAOVA', 'DSLDANSLT', 'medLDA', 'DSLDA', 'DSLDAOSST', 'random'}, 'Location','NorthWest'); 
    hold off;
    saveas(gcf,[savefilepath num2str(minvtopic) '_' titlestr '.jpg']);
end

end

