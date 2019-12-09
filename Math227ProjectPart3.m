npoints=200;
me=9.10938188e-31; %mass of electron 
length=10e-9;   %10nm
x=0:length/npoints:length;

hbar=1.054571596;										%Planck's constant (x10^34 J s)

n=4  % up to first 4 bound states
s1=char('.y','.k','.r','.g','.b','.m','.c');
s=char('b','r','y','m','b','m','c');								%plot curves in different colors


%case 1 no perubation________________________________________________



%potential well diagram
figure();
for i=1:npoints+1; v(i)=0; end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 ttl=['CASE 1, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 title(ttl);

%numerical solution 

[energy,phi]=solve_schM(10,npoints,v,1,n);		%call solve_schM
figure()
for i=1:4
       
    
    %j=1+mod(i,7);
      plot(x,(phi(:,i)).^2,s1(i));											%plot eigenfunctions
       hold on;
       sprintf(['eigenenergy CASE 1 (',num2str(i),') = ',num2str(energy(i)),' eV'])		%energy eigenvalues
end
 tt3=['CASE 1, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Probability density');
 title(tt3);
 
 %CASE 2:pertubation_____________________________
 
 %potential well diagram
figure();
for i=1:npoints+1
    if i<npoints/2+1
        v(i)=0;
    else
        v(i)=0.01;
    end
    
     end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 tt4=['CASE 2, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 title(tt4);
 
 %numerical solution 

[energy,phi]=solve_schM(10,npoints,v,1,n);		%call solve_schM
figure()
for i=1:4
       
    
    %j=1+mod(i,7);
      plot(x,(phi(:,i)).^2,s1(i));											%plot eigenfunctions
       hold on;
       sprintf(['eigenenergy CASE 2(',num2str(i),') = ',num2str(energy(i)),' eV'])		%energy eigenvalues
end
 tt4=['CASE 2,Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Probability density');
 title(tt4);
 
  %CASE 3____________________________________________________

%potential well diagram
figure();
for i=1:npoints+1
    if i<(npoints/4+1) | ((i<3*npoints/4+1) &(i>npoints/2+1)) 
        v(i)=0.01;
    else
        v(i)=0.0;

    end
    
     end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 tt4=['CASE 3, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 title(tt4);
 
 %numerical solution 

[energy,phi]=solve_schM(10,npoints,v,1,n);		%call solve_schM
figure()
for i=1:4
       
    
    %j=1+mod(i,7);
      plot(x,(phi(:,i)).^2,s1(i));											%plot eigenfunctions
       hold on;
       sprintf(['eigenenergy CASE 3(',num2str(i),') = ',num2str(energy(i)),' eV'])		%energy eigenvalues
end
 tt4=['Potential Well(CASE 3),Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Probability density');
 title(tt4);
 
 
 %CASE 4____________________________________________________

%potential well diagram
figure();
for i=1:npoints+1
    if i<(npoints/4+1) | ((i>3*npoints/4+1)) 
        v(i)=0.01;
    else
        v(i)=0.0;

    end
    
     end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 tt4=['CASE 4, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 title(tt4);
 
 %numerical solution 

[energy,phi]=solve_schM(10,npoints,v,1,n);		%call solve_schM
figure()
for i=1:4
       
    
    %j=1+mod(i,7);
      plot(x,(phi(:,i)).^2,s1(i));											%plot eigenfunctions
       hold on;
       sprintf(['eigenenergy CASE 4(',num2str(i),') = ',num2str(energy(i)),' eV'])		%energy eigenvalues
end
 tt4=['CASE 4,Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Probability density');
 title(tt4);
 
 %CASE 5____________________________________________________

%potential well diagram
figure();
for i=1:npoints+1
    if i>(npoints/4+1) & (i<3*npoints/4+1) 
        v(i)=0.01;
    else
        v(i)=0.0;

    end
    
     end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 tt4=['CASE 5, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 title(tt4);
 
 %numerical solution 

[energy,phi]=solve_schM(10,npoints,v,1,n);		%call solve_schM
figure()
for i=1:4
       
    
    %j=1+mod(i,7);
      plot(x,(phi(:,i)).^2,s1(i));											%plot eigenfunctions
       hold on;
       sprintf(['eigenenergy CASE 5(',num2str(i),') = ',num2str(energy(i)),' eV'])		%energy eigenvalues
end
 tt4=['CASE 5, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Probability density');
 title(tt4);
 
 %CASE 6____________________________________________________

%potential well diagram
figure();
for i=1:npoints+1
        v(i)=0.01*sin(i);
    end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 tt4=['CASE 6, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 title(tt4);
 
 %numerical solution 

[energy,phi]=solve_schM(10,npoints,v,1,n);		%call solve_schM
figure()
for i=1:4
       
    
    %j=1+mod(i,7);
      plot(x,(phi(:,i)).^2,s1(i));											%plot eigenfunctions
       hold on;
       sprintf(['eigenenergy CASE 6(',num2str(i),') = ',num2str(energy(i)),' eV'])		%energy eigenvalues
end
 tt4=['CASE 6,Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Probability density');
 title(tt4);
 
%  (2*(((2*pi*a*l*n-a^2*l^2)*cos(2*pi*n+a*l+1)+(-2*pi*a*l*n-a^2*l^2)*cos(2*pi*n-a*l-1)-8*pi^2*cos(a*l+1)*n^2+2*a^2*l^2*cos(a*l+1))/(16*pi^2*a*n^2-4*a^3*l^2)+(2*pi^2*cos(1)*n^2)/(4*pi^2*a*n^2-a^3*l^2))*v)/l
 
 
  %CASE 7____________________________________________________

%potential well diagram
figure();
for i=1:npoints+1
        v(i)=0.01*sin(i/20);
end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 tt4=['CASE 7, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 title(tt4);
 
 %numerical solution 

[energy,phi]=solve_schM(10,npoints,v,1,n);		%call solve_schM
figure()
for i=1:4
       
    
    %j=1+mod(i,7);
      plot(x,(phi(:,i)).^2,s1(i));											%plot eigenfunctions
       hold on;
       sprintf(['eigenenergy CASE 7(',num2str(i),') = ',num2str(energy(i)),' eV'])		%energy eigenvalues
end
 tt4=['CASE 7, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Probability density');
 title(tt4);
 
 %(2*(((200*pi*a*l*n-5*a^2*l^2)*cos((40*pi*n+a*l+1)/20)+(-200*pi*a*l*n-5*a^2*l^2)*cos((40*pi*n-a*l-1)/20)-16000*pi^2*cos((a*l+1)/20)*n^2+10*a^2*l^2*cos((a*l+1)/20))/(1600*pi^2*a*n^2-a^3*l^2)+(16000*pi^2*cos(1/20)*n^2)/(1600*pi^2*a*n^2-a^3*l^2))*v)/l
 
  %CASE 8____________________________________________________

%potential well diagram
figure();
for i=1:npoints+1
        v(i)=0.01*cos(i/20);
    end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 tt4=['CASE 8, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 title(tt4);
 
 %numerical solution 

[energy,phi]=solve_schM(10,npoints,v,1,n);		%call solve_schM
figure()
for i=1:4
       
    
    %j=1+mod(i,7);
      plot(x,(phi(:,i)).^2,s1(i));											%plot eigenfunctions
       hold on;
       sprintf(['eigenenergy CASE 8(',num2str(i),') = ',num2str(energy(i)),' eV'])		%energy eigenvalues
end
 tt4=['CASE 8, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Probability density');
 title(tt4);
 
 %(2*(-((200*pi*a*l*n-5*a^2*l^2)*sin((40*pi*n+a*l+1)/20)+(200*pi*a*l*n+5*a^2*l^2)*sin((40*pi*n-a*l-1)/20)-16000*pi^2*sin((a*l+1)/20)*n^2+10*a^2*l^2*sin((a*l+1)/20))/(1600*pi^2*a*n^2-a^3*l^2)-(16000*pi^2*sin(1/20)*n^2)/(1600*pi^2*a*n^2-a^3*l^2))*v)/l
 
   %CASE 9____________________________________________________

%potential well diagram
figure();
for i=1:npoints+1
    if i<npoints/2
        v(i)=0.1-0.01*i/(10);
    else
        v(i)=-0.1+0.01*i/(10);
    end
        
    end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 tt4=['CASE 9, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 title(tt4);
 
 %numerical solution 

[energy,phi]=solve_schM(10,npoints,v,1,n);		%call solve_schM
figure()
for i=1:4
       
    
    %j=1+mod(i,7);
      plot(x,(phi(:,i)).^2,s1(i));											%plot eigenfunctions
       hold on;
       sprintf(['eigenenergy CASE 9(',num2str(i),') = ',num2str(energy(i)),' eV'])		%energy eigenvalues
end
 tt4=['CASE 9, m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Probability density');
 title(tt4);





