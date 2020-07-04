% Creates the plots for the equally weighted portfolios.
function c = logcumplot(returns, titles)
a = log(cumprod(returns/100 + 1));
b = csvread("12_indust_month_value.csv",1 ,0);
b2 = b(end-length(returns)+1:end,1);
c = a(end);
time = datenum(num2str(b2),'yyyymm');
annRet = num2str(round(mean(returns)*12,2));
annRetStd = strcat("[", num2str(round(sqrt(12)*std(returns),2)),"]");
MaxDD = num2str(round(maxdrawdown(cumprod(returns/100 + 1))*100,2));
figure(gcf);
plot(time,a,'b', 'LineWidth', 2);hold on;
plot([min(xlim()),max(xlim())],[0,0], 'k--');
ax = gca; 
txt1 = strcat("Ann. average return = ",annRet,"% ", annRetStd);
txt2 = strcat("Max Drawdown =  ",MaxDD, "%");
text(time(5),3.30,txt1,'FontSize',13)
text(time(5),2.95,txt2,'FontSize',13)
ylim([-0.5, 3.5]);
xticks(time(1:85:end));
datetick('x', 'yyyy', 'keepticks')
ax.FontSize = 12;
recessionplot;
set(gcf,'units','points','position',[10,10,400,250])
end