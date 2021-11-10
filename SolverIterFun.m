clear all

h=0.1;
dw=0.5;

ff=[];

Ldw=1.5;
dSOC=1;

aaa=15:15

signSOC=1;

h1=figure
h2=figure

for indN=2:2
    
    indN
    if (indN==1) 
        N1=90;
  indT=64;
    end 
    
    if (indN==2)   
        N1=50;
 indT=68;
    end 
    
    if (indN==3) 
      N1=20;
 indT=115;
    end 
    
    fileName=['FunSO' num2str(aaa) 'dSoc' num2str(dSOC) 'indN' num2str(indN)];
            
N=2*N1;

IndJ=1;
IndJPrev=IndJ;

ind=0;

%for indT=21:5:150
    ind=ind+1;
indT
beta=10;
am=0.004;
ep=-1;

tau=0.0025*(indT-0.9)+0.001;
   
J=-am*beta;
a=beta*(Ldw)^2*(-tau); 
b=beta*(Ldw)^2*0.023;
c=0.37*beta; 
%% change sign of alphaSO
alphaSO= -signSOC*2*(0.1*aaa*Ldw)*am*beta;

Kan=-1.0*beta*(0.1*aaa*Ldw)^2*am/dSOC;

x=h*(1:N) - h*N1;
L=h*(N) - h*N1;

%if(ind==1)
theta0= pi* ( tanh(x/dw)/abs(tanh(L/dw))  +0)/2 ;
psi0=(-(a+Kan)/(2*b))^(0.5)  +0*x;
%end

dev=1;

(1-J*a/(2*b));
Flag=1;
while((dev>0.0001) )
      
psiPrev=psi0;

[theta] = FunSolveLLG(h,N,theta0,psi0,a,b,c,alphaSO,J,Kan,ep);
theta0=theta;

%plot(x,theta);

%pause
psi = FunSolveGL(h,N,theta0,psi0,a,b,c,alphaSO,J,Kan);
psi0=psi;

%plot(x,psi)
dev=max(abs(psiPrev-psi0))

%pause
    
end

psi=(psi(1:length(x)) + psi(length(x):-1:1))/2;

%%%%%%%%%%%%%%%%%%%%%%% calculate energy %%%%%%%%%%%

figure(h1)

plot(x/max(x),psi)
hold on

figure(h2)
plot(x/max(x),theta)
hold on

end

save(fileName,'x', 'psi', 'theta');
%save(fileName,'ttau', 'Eenergy1', 'Eenergy2','EenergyDW');
%load(fileName,'ttau', 'NNN', 'Eenergy');

%L=(NNN*h);
% figure
% 
% plot(ttau,EenergyDW-Eenergy1)
% hold on
% plot(ttau,EenergyDW-Eenergy2)

%plot(L,Eenergy./L)
%  %[xq,yq] = meshgrid( min(L):h/2:max(L), min(ttau):0.01/2:max(ttau));
%  %vq = griddata(L,ttau,Eenergy,xq,yq);
% % 
%  % surf(xq,yq,vq./xq)
%  %   shading interp
%    
% %    %%%%%%%%%%%%%%%%%%%%%%%%
% %    
% %    figure
% % %    
%      [xq,yq] = meshgrid( min(L):h/2:max(L), -min(ttau):-0.01/2:-max(ttau));
%   vq = griddata(L,-ttau,Eenergy,xq,yq);
%      [hh,hh]=contourf(xq,yq,vq./xq,1000)
% % %   % [hh,hh]=contourf(xq,yq,vq,100)
%    set(hh,'edgecolor','none');
% colormap(bluewhitered(1000))
% %  %  [hh,hh]=contourf(xq,yq,vq./xq,1000)
%   colorbar
% %  
% %   ylabel('$T/T_c-1$','interpreter','latex','FontSize',26)
% %  xlabel('$L_{dw}/l_{Bl}$','interpreter','latex','FontSize',26)
% %  
% %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  
%  TTmin=[];
%  LLmin=[];
%  
%  TT0=[];
%  LL0=[];
%  vvq=vq./xq;
%  for ii=1:length(xq(:,1))
% 
%              [E0,I]= min(vvq(ii,:).^2);
%              
%              if( yq(ii,I)< -0.02 )
%     
% LL0=[LL0 xq(ii,I)];
%         TT0=[TT0 yq(ii,I)];
%       end
%         
%         [Emin,I]= min(vvq(ii,:));
%     if(Emin<0)
%     
%         LLmin=[LLmin xq(ii,I)];
%         TTmin=[TTmin yq(ii,I)];
%        
%       end
%     
%  end
%  
%  
%   plot(LLmin,TTmin, 'linewidth',2,'color',[1,0.65,0])
%   hold 
%   plot(LL0,TT0,'linewidth',2,'color','k')
% %   
% % fname1=[fileName '.png']
% %   print(gcf,fname1,'-dpng','-r300')



