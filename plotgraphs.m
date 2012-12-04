function [] = plotgraphs()
% UNTITLED2 Summary of this function goes here
% Detailed explanation goes here

minvtopic = 1;
filename1 = ['finalbigrun6.mat'];
filename2 = ['finalbigrun6_std.mat'];
savefilepath = '/v/filer4b/v16q001/data/Ayans/Results_finalbigrun6/';

A = load(filename1);
B = load(filename2);
k2mat = [20:10:200];

for k2ind=1:length(k2mat)
    
    C = squeeze(A.savedstat(:,k2ind,:));
    C(:,3) = C(:,3); 
    D = squeeze(B.savedstatstd(:,k2ind,:));
    p = [0.1:0.1:1]';
    rbaseline = 0.23*ones(size(p));
    
    figure (1);
    
    %errorbar(p,C(:,1),D(:,1),'r.-','LineWidth',2, 'MarkerEdgeColor','r','MarkerSize',15); hold on;
    errorbar(p,C(:,2),D(:,2),'g.-','LineWidth',2, 'MarkerEdgeColor','g','MarkerSize',15); hold on;
    errorbar(p,C(:,3),D(:,3),'b.-','LineWidth',2, 'MarkerEdgeColor','b','MarkerSize',15); hold on;
    errorbar(p,C(:,4),D(:,4),'k.-','LineWidth',2, 'MarkerEdgeColor','k','MarkerSize',15); hold on;
    errorbar(p,C(:,5),D(:,5),'r.-','LineWidth',2, 'MarkerEdgeColor','r','MarkerSize',15); hold on;
    errorbar(p,C(:,6),D(:,6),'m.-','LineWidth',2, 'MarkerEdgeColor','m','MarkerSize',15); hold on;
    plot(p,rbaseline,'c.-','LineWidth',2); grid on, hold on;
    t = xlabel('fraction of training data'); hold on;
    set(t, 'FontSize', 16);
    t = ylabel('over-all classification accuracy'); hold on;
    set(t, 'FontSize', 16);
    titlestr = ['Learning curve for k2 = ' num2str(k2mat(k2ind))];
    t = title(titlestr);
    set(t, 'FontSize', 16);
    set(findobj('type','axes'),'fontsize',14);
    set(gca,'YLim',[0.1 1.0]);
    set(gca,'XLim',[0 1.1]);
    hold on;
    
    legend({'medLDA-OVA', 'DSLDA-NSLT', 'medLDA-MTL', 'DSLDA', 'DSLDA-OSST', 'random'}, 'Location','NorthWest'); 
    %%legend boxoff; 
    hold off;
    saveas(gcf,[savefilepath num2str(minvtopic) '_' titlestr '.jpg']);
end

end

