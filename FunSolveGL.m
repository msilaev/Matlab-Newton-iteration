% Equation d^2_x \psi = a \psi + b \psi^3

%clear all

% J=0.04*(10/Ldw)^2;
% a=10^3*(0.026-tau); 
% a=10^3*(-tau); 
% 
% b=10^3*0.023;
% c=0.04*10^3*(10/Ldw)^2; 
% alphaSO=-10*0.04*20*0.1*(10/Ldw);


function psi = FunSolveGL(h,N,theta,psi0,a,b,c,alphaSO,J,Kan)

A = zeros(N,N); 
Rhs=zeros(N,1);

psi=psi0;
  for i=2:N-1
   
    A(i,i) = 4*c/h^2 + 2*a + 12*b*psi(i)^2 + 2*J*( theta(i) - theta(i-1) )^2/h^2 + 2*alphaSO*(theta(i) - theta(i-1))/h +...
    2*Kan*(cos(theta(i)))^2;
    A(i,i-1) = -2*c/h^2; 
    A(i,i+1) = -2*c/h^2;
    
    Rhs(i) = 2*c*(2*psi(i) - psi(i-1) - psi(i+1) )/h^2 + 2*a*psi(i) + 4*b*psi(i)^3 + ...
        ( 2*J* (theta(i) - theta(i-1) )^2/h^2 + 2*alphaSO* ( theta(i) - theta(i-1) )/h  +  2*Kan*(cos(theta(i)))^2 )*psi(i) ;
end

A(1,1)=(2*c/h^2 + 2*a + 12*b*psi(1)^2) + 2*J*( theta(1) +pi/2 )^2/h^2 + 2*alphaSO*(theta(1) +pi/2)/h + 2*Kan*(cos(theta(1)))^2; 
A(1,2)=-2*c/h^2;

A(N,N)=(2*c/h^2 + 2*a + 12*b*psi(N)^2)+ 2*J*( pi/2 - theta(N) )^2/h^2 + 2*alphaSO*(pi/2-theta(N))/h + 2*Kan*(cos(theta(N)))^2; 
A(N-1,N)=-2*c/h^2;

Rhs(1) = 2*c*(psi(1) - psi(2) )/h^2 + 2*a*psi(1) + 4*b*psi(1)^3 + ...
        ( 2*J* (theta(1) +pi/2 )^2/h^2 + 2*alphaSO* ( theta(1) + pi/2 )/h  + 2*Kan*(cos(theta(1)))^2 )*psi(1) ;
Rhs(N) = 2*c*(psi(N) - psi(N-1) )/h^2 + 2*a*psi(N) + 4*b*psi(N)^3 + ...
        ( 2*J* (pi/2-theta(N)  )^2/h^2 + 2*alphaSO* ( pi/2-theta(N)  )/h + 2*Kan*(cos(theta(N)))^2 )*psi(N) ;


while ( max( abs(Rhs) )> 0.0001)
        
    max( abs(Rhs) );
         
    psi = (psi' - inv(A)*Rhs)';
    %psi(N)=psi(1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   for i=2:N-1
   
    A(i,i) = 4*c/h^2 + 2*a + 12*b*psi(i)^2 + 2*J*( theta(i) - theta(i-1) )^2/h^2 + 2*alphaSO*(theta(i) - theta(i-1))/h + ...
    2*Kan*(cos(theta(i)))^2; 

    A(i,i-1) = -2*c/h^2; 
    A(i,i+1) = -2*c/h^2;
    
    Rhs(i) = 2*c*(2*psi(i) - psi(i-1) - psi(i+1) )/h^2 + 2*a*psi(i) + 4*b*psi(i)^3 + ...
        ( 2*J* (theta(i) - theta(i-1) )^2/h^2 + 2*alphaSO* ( theta(i) - theta(i-1) )/h + 2*Kan*(cos(theta(i)))^2 )*psi(i) ;
end

A(1,1)=(2*c/h^2 + 2*a + 12*b*psi(1)^2) + 2*J*( theta(1) +pi/2 )^2/h^2 + 2*alphaSO*(theta(1) + pi/2)/h + 2*Kan*(cos(theta(1)))^2;  
A(1,2)=-2*c/h^2;

A(N,N)=(2*c/h^2 + 2*a + 12*b*psi(N)^2)+ 2*J*(theta(N) - theta(N-1) )^2/h^2 + 2*alphaSO*(theta(N) - theta(N-1))/h + ...
2*Kan*(cos(theta(N)))^2;
A(N-1,N)=-2*c/h^2;

Rhs(1) = 2*c*(psi(1) - psi(2) )/h^2 + 2*a*psi(1) + 4*b*psi(1)^3 + ...
        ( 2*J* (theta(1) +pi/2 )^2/h^2 + 2*alphaSO* ( theta(1) +pi/2 )/h + 2*Kan*(cos(theta(1)))^2 )*psi(1)  ;
    
Rhs(N) = 2*c*(psi(N) - psi(N-1) )/h^2 + 2*a*psi(N) + 4*b*psi(N)^3 + ...
        ( 2*J* (theta(N)-theta(N-1)  )^2/h^2 + 2*alphaSO* ( theta(N)-theta(N-1) )/h + 2*Kan*(cos(theta(N)))^2 )*psi(N) ;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


end



%h1=figure;
%figure(h1)

%plot(x,psiPrev)
 %   pause