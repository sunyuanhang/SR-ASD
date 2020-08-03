close all
clc
clear
load fz1.mat %x3 is fault data.
load fz1_before_adding_the_noise.mat
NN=2048;% the length of the signal
fs=10000; % sampling rate
t=1/fs;%sampling time
%% original simulation data
figure(1)
plot(x1);
axis([0,2048,-inf,inf]);
set(gcf,'Position',[100 100 380 190]);
set(gca,'Position',[.13 .17 .80 .74]);
figure_FontSize=10;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',1);
str=strcat('Amplitude/(m/s^2)');
xlabel('Sample Number','FontSize',11,'Fontname','Times New Roman');
ylabel(str,'FontSize',11,'Fontname','Times New Roman')
%%
figure(2)
plot(x3);
axis([0,2048,-inf,inf]);
set(gcf,'Position',[100 100 380 190]);
set(gca,'Position',[.13 .17 .80 .74]);
figure_FontSize=10;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',1);
str=strcat('Amplitude/(m/s^2)');
xlabel('Sample Number','FontSize',11,'Fontname','Times New Roman');
ylabel(str,'FontSize',11,'Fontname','Times New Roman')

%% SR-ASD method
y=x3';
lam1=3.3514
lam2 =1.1282  % lam1 and lam2 are set based on the noise standard deviation
[x]=SR_ASD(y,lam1,lam2);
figure (4)
plot(x);
axis([0,2048,-inf,inf]);
set(gcf,'Position',[100 100 380 190]);
set(gca,'Position',[.13 .17 .80 .74]);
figure_FontSize=10;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',1);
str=strcat('Amplitude/(m/s^2)');
 xlabel('Sample Number','FontSize',11,'Fontname','Times New Roman');
 ylabel(str,'FontSize',11,'Fontname','Times New Roman')
%% FFT
ss=x;
yyy=abs(ss);
 y_hht=hilbert(ss);
y_an=abs(y_hht);
y_an=y_an-mean(y_an);
y_an_nfft= 2^nextpow2(length(y_an));%
y_an_ft=fft(y_an,y_an_nfft);%
y_an_f=fs*(0:y_an_nfft/2-1)/y_an_nfft;%
y_an_p=y_an_ft.*conj(y_an_ft)/y_an_nfft;
c=2*abs(y_an_ft(1:y_an_nfft/2))/length(y_an);
d=y_an_p(1:y_an_nfft/2);
N=length(c);
figure(5)
plot(y_an_f(1:N),c(1:N));    
set(gcf,'Position',[100 100 380 190]);
set(gca,'Position',[.16 .20 .80 .74]);
figure_FontSize=10;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',1);
 str=strcat('Amplitude/(m/s^2)');
 xlabel('Frequency/Hz','FontSize',10,'Fontname','Times New Roman');
 ylabel(str,'FontSize',10,'Fontname','Times New Roman')

 
