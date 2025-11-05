%function f_Main_Figure
% figure 5(c),(d) 
close all
k=7;t=2;
%k= 6

Res= log(1./Response(:,k));
KACF2= KurtosisStdR1(t,:)';
Rec= log(1./Recovery(:,k));

X= KACF2;
Y= Res;
[CC,PP]=corr(X,Y,'Type','Spearman')
P = polyfit(X,Y,1);
yfit = P(1)*X+P(2);

GDP1= rescale(GDP,100,1000);
SZ= floor(GDP1);
SZe= SZ(Eindx);SZd= SZ(Dindx);

figure;
h3= scatter( KACF2(Eindx),Res(Eindx),SZe,'r','filled','MarkerFaceAlpha',0.7, 'MarkerEdgeColor', [0.8, 0.8, 0.8],'linewidth',2.5); hold on;
h4= scatter(KACF2(Dindx), Res(Dindx),SZd,'b','filled','MarkerFaceAlpha',0.7, 'MarkerEdgeColor', [0.8, 0.8, 0.8],'linewidth',2.5); 
[Value, Index] = sort(X);plot(X(Index),yfit(Index),'k--', 'LineWidth',3, 'Color', [0.4, 0.4, 0.4]);
h1_legend = plot(nan, nan, 'ro', 'MarkerSize', 20, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', [0.8, 0.8, 0.8],'linewidth',2.5');
h2_legend = plot(nan, nan, 'bo', 'MarkerSize', 20, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', [0.8, 0.8, 0.8],'linewidth',2.5);
legend([h1_legend, h2_legend], 'Emerging countries', 'Developed countries', 'Box', 'off');
xlabel('KACF','fontsize',14);
ylabel('Response time','fontsize',14);
set(gca,'fontsize',18);
set(gcf, 'Color', 'w');
set(gca,'box','on', 'linewidth',3)
print('Fig5c_ResposneTime_Window120_alpha100.png', '-dpng', '-r1200'); 


Y= Rec;
P = polyfit(X,Y,1);
yfit = P(1)*X+P(2);
[CC,PP]=corr(X,Y,'Type','Spearman')

figure;
h3= scatter( KACF2(Eindx),Rec(Eindx),SZe,'r','filled','MarkerFaceAlpha',0.7, 'MarkerEdgeColor', [0.8, 0.8, 0.8],'linewidth',2.5); hold on;
h4= scatter(KACF2(Dindx), Rec(Dindx),SZd,'b','filled','MarkerFaceAlpha',0.7, 'MarkerEdgeColor', [0.8, 0.8, 0.8],'linewidth',2.5); 
[Value, Index] = sort(X);plot(X(Index),yfit(Index),'k--', 'LineWidth',3, 'Color', [0.4, 0.4, 0.4]);
xlabel('KACF','fontsize',14); ylabel('Recovery time','fontsize',14);
set(gca,'fontsize',18); set(gcf, 'Color', 'w'); set(gca,'box','on', 'linewidth',3)

print('Fig5d_ResposneTime_Window120_alpha100.png', '-dpng', '-r1200'); 