% Equation d^2_x \psi = a \psi + b \psi^3

%clear all

% J=0.04*(10/Ldw)^2;
% a=10^3*(0.026-tau); 
% a=10^3*(-tau); 
% 
% b=10^3*0.023;
% c=0.04*10^3*(10/Ldw)^2; 
% alphaSO=-10*0.04*20*0.1*(10/Ldw);


function [theta] = FunSolveLLG(h,N,theta0,psi,a,b,c,alphaSO,J,Kan,ep)

A = zeros(N,N); 
Rhs=zeros(N,1);

theta=theta0;
for i=2:N-1
   
    A(i,i) = 4/h^2 + 2*(ep-Kan*psi(i)^2)*cos(2*theta(i)) + 2*J*( psi(i)^2 + psi(i+1)^2 )/h^2 ;
    A(i,i-1) = -2*(1+J*psi(i)^2)/h^2; 
    A(i,i+1) = - 2*(1+J*psi(i+1)^2)/h^2;
    
    Rhs(i) = 2*(2 + J*(psi(i)^2 + psi(i+1)^2))*theta(i)/h^2 ...
        - 2*(1 + J*psi(i)^2)*theta(i-1)/h^2 ...
        - 2*(1 + J*psi(i+1)^2)*theta(i+1)/h^2 ...
        + (ep-Kan*psi(i)^2)*sin(2*theta(i)) ...
        + alphaSO* ( psi(i)^2 - psi(i+1)^2 )/h ;
    
end
    
A(1,1)=2*(2 + J*( psi(1)^2 + psi(2)^2 ))/h^2 +2*(ep-Kan*psi(1)^2)*cos(2*theta(1)) ;
A(1,2)=-2*(1+J*psi(2)^2)/h^2;

A(N,N)=4*(1 + J*psi(N)^2 )/h^2 + 2*(ep-Kan*psi(N)^2)*cos(2*theta(N)) ;
A(N-1,N)=-2*(1+J*psi(N)^2)/h^2;

% Rhs(1) = 2*(2 + J*(psi(1)^2 + psi(2)^2))*(theta(1) + pi/2)/h^2 ...
%         - 2*(1 + J*psi(2)^2)*theta(2)/h^2 ...
%         + (ep-Kan*psi(1)^2)*sin(2*theta(1)) ...
%         + alphaSO* ( psi(1)^2 - psi(2)^2 )/h ;

Rhs(1) = 2*(1 + J*(psi(1)^2 + psi(2)^2)/2)*(2*theta(1)- theta(2) + 1*pi/2)/h^2 ...
                + (ep-Kan*psi(1)^2)*sin(2*theta(1)) ...
        + alphaSO* ( psi(1)^2 - psi(2)^2 )/h ;

Rhs(N) = 2*(1 + J*(psi(N)^2))*(2*theta(N) - theta(N-1) - 1*pi/2)/h^2 ...
        + (ep-Kan*psi(N)^2)*sin(2*theta(N)) ;
        
%     Flag=0;
%     if( abs(det(A))>1 )
%         Flag=1 ;
%     end
    
while ( max( abs(Rhs) )> 0.0001 )
        
    max( abs(Rhs) );
    
             
   % theta = (theta' - inv(A)*Rhs)';
    theta =0.999*theta + 0.001* (theta' - inv(A)*Rhs)';
    
%     plot(theta)
%      pause
%     hold on
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   for i=2:N-1
   
    A(i,i) = 4/h^2 + 2*(ep-Kan*psi(i)^2)*cos(2*theta(i)) + 2*J*( psi(i)^2 + psi(i+1)^2 )/h^2 ;
    A(i,i-1) = -2*(1+J*psi(i)^2)/h^2; 
    A(i,i+1) = - 2*(1+J*psi(i+1)^2)/h^2;
    
    Rhs(i) = 2*(2 + J*(psi(i)^2 + psi(i+1)^2))*theta(i)/h^2 ...
        - 2*(1 + J*psi(i)^2)*theta(i-1)/h^2 ...
        - 2*(1 + J*psi(i+1)^2)*theta(i+1)/h^2 ...
        + (ep-Kan*psi(i)^2)*sin(2*theta(i)) ...
        + alphaSO* ( psi(i)^2 - psi(i+1)^2 )/h ;
    
end
    
A(1,1)=2*(2 + J*( psi(1)^2 + psi(2)^2 ))/h^2 +2*(ep-Kan*psi(1)^2)*cos(2*theta(1)) ;
A(1,2)=-2*(1+J*psi(2)^2)/h^2;

A(N,N)=4*(1 + J*psi(N)^2 )/h^2 + 2*(ep-Kan*psi(N)^2)*cos(2*theta(N)) ;
A(N-1,N)=-2*(1+J*psi(N)^2)/h^2;

Rhs(1) = 2*(1 + J*(psi(1)^2 + psi(2)^2)/2)*(2*theta(1)-theta(2)+1*pi/2)/h^2 ...
        - 0*(1 + J*psi(2)^2)*theta(2)/h^2 ...
        + (ep-Kan*psi(1)^2)*sin(2*theta(1)) ...
        + alphaSO* ( psi(1)^2 - psi(2)^2 )/h ;

Rhs(N) = 2*(1 + J*(psi(N)^2))*(2*theta(N) - theta(N-1) -1*pi/2)/h^2 ...
        + (ep-Kan*psi(N)^2)*sin(2*theta(N)) ;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%     if( abs(det(A))<1 )
%         Flag=0 ;
%     end

end


    end



    