clear all
close all
clc

%get B6 from batch 0

filename = []

b = importdata('filtered_sorted_otu_table_gg_w_taxon_L6.txt');
exp1_data_ids = b.textdata(1,2:end)';
exp1_data = b.data;
exp1_data = exp1_data./repmat(sum(exp1_data),size(exp1_data,1),1);
exp1_bugs = b.textdata(2:end,1);

exp1_hc_data_index = startsWith(exp1_data_ids,'HC.F','IgnoreCase',true) | startsWith(exp1_data_ids,'HC.KO.F','IgnoreCase',true)
exp1_glut_data_index = startsWith(exp1_data_ids,'Glut.F','IgnoreCase',true) | startsWith(exp1_data_ids,'Glut.KO.F','IgnoreCase',true)
exp1_chow_data_index = startsWith(exp1_data_ids,'Chow.F','IgnoreCase',true) | startsWith(exp1_data_ids,'Chow.KO.F','IgnoreCase',true)

bit = exp1_hc_data_index | exp1_glut_data_index | exp1_chow_data_index;

[coeff,SCORE,latent,tsquared,explained]  = pca((exp1_data(:,bit)'))
colors = [1 0 0; 0 0 0; 0 0 1];

x1 = startsWith(exp1_data_ids(bit),'HC.F','IgnoreCase',true) |startsWith(exp1_data_ids(bit),'HC.KO.F','IgnoreCase',true) 
x2 = startsWith(exp1_data_ids(bit),'Glut.F','IgnoreCase',true) | startsWith(exp1_data_ids(bit),'Glut.KO.F','IgnoreCase',true)
x3 = startsWith(exp1_data_ids(bit),'Chow.F','IgnoreCase',true) | startsWith(exp1_data_ids(bit),'Chow.KO.F','IgnoreCase',true)

%c1 = lines(4);
%colors = [c1(3,:);c1(1,:);c1(2,:);c1(4,:)];


figure
scatter(SCORE(x1,1),SCORE(x1,2),150,repmat(colors(1,:),sum(x1),1),'o','filled')
hold on
pause
scatter(SCORE(x2,1),SCORE(x2,2),150,repmat(colors(2,:),sum(x2),1),'o','filled')
pause
scatter(SCORE(x3,1),SCORE(x3,2),150,repmat(colors(3,:),sum(x3),1),'o','filled')


grid off; %// Show a grid
set(0,'defaulttextinterpreter','none')
set(gca,'Box','off','FontName','Helvetica','FontSize', 25,'fontweight','bold','LineWidth',2)
ylabel(''), xlabel(''), zlabel(''),title('');

explained(1)
explained(2)

pause


figure
sum((exp1_data(:,bit)'))
[q,w] = sort(mean((exp1_data(:,bit)')))
tmp = [exp1_data(w(end-4:end),bit); ...
      1-sum(exp1_data(w(end-4:end),bit))]
  
tmp2 = [tmp(:,x1), tmp(:,x2), tmp(:,x3)]  
  tt = exp1_bugs(w(end-4:end));
  tt(end:-1:1)

bar(tmp2', 'stacked','FaceColor','flat')
grid off
set(0,'defaulttextinterpreter','none')
set(gca,'XGrid','off','YGrid','off','Box','off','FontName','Helvetica','FontSize', 25,'fontweight','bold','LineWidth',2)
ylabel(''), xlabel(''), title('');



