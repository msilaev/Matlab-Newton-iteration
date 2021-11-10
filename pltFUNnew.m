clear all

Ldw=1.5;
dSOC=1;

aaa=15:15

signSOC=1;

h1=figure
h2=figure

for indN=1:3
    
    if (indN==1) 
        N1=90;
  indT=64;
  ccolor='blue';
    end 
    
    if (indN==2)   
        N1=50;
 indT=68;
 
 ccolor='green';
    end 
    
    if (indN==3) 
      N1=20;
 indT=115;
 ccolor='red';
    end 
            
    tau=0.0025*(indT-0.9)+0.001;
    
a=10^3*(-tau); 

b=10^3*0.023;
    
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    fileName=['FunSO' num2str(aaa) 'dSoc' num2str(dSOC) 'indN' num2str(indN)];
    load(fileName,'x', 'psi', 'theta');
    
    psi0=(-a/(2*b))^(0.5)  +0*x;
    
figure(h1)
plot(x/max(x),psi./psi0-1,'LineWidth',2,'color',ccolor)
hold on

figure(h2)
plot(x/max(x),theta/pi,'LineWidth',2,'color',ccolor)
hold on

end

figure(h1)
set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',26)
  
ylabel('$\psi/\psi_0-1$','interpreter','latex','FontSize',28)
xlabel('$x/L_{pdw}$','interpreter','latex','FontSize',28)

grid on

fileNameFig=[fileName 'FigPsi'];

fname1=[fileNameFig '.png']
%set(gcf, 'PaperUnits', 'inches');
%set(gcf, 'PaperSize', [10 10]);
print(gcf,fname1,'-dpng','-r300')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(h2)
set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',26)
  
ylabel('$\theta/\pi$','interpreter','latex','FontSize',28)
xlabel('$x/L_{pdw}$','interpreter','latex','FontSize',28)

grid on

fileNameFig=[fileName 'FigTheta'];
fname1=[fileNameFig '.png']
%set(gcf, 'PaperUnits', 'inches');
%set(gcf, 'PaperSize', [10 10]);
print(gcf,fname1,'-dpng','-r300')



save(fileName,'x', 'psi', 'theta');
