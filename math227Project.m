%Infinite potential well

%Compare QM (time-independant)pertubation theory to approximate 
%the Energy levels and wave function from Schrodinger's eq

%compare with mathematical techniques

%consider a potential well with small pertubation
%hamiltonian H=H0+lambda*Hprime

%En=<siN|Hprime|siN>    (1st order)
%siN=sum(n!=m) (<sim0|H0|sin0>/(En0-Em0))siM0   
clear
clf;
%CASE 1:no pertubation___________________

%analytical solution

npoints=200;
me=9.10938188e-31; %mass of electron 
length=10e-9;   %10nm
x=0:length/npoints:length;

hbar=1.054571596;										%Planck's constant (x10^34 J s)

n=4  % up to first 4 bound states
s1=char('.y','.k','.r','.g','.b','.m','.c');
s=char('b','r','y','m','b','m','c');								%plot curves in different colors

%potential well diagram
figure();
for i=1:npoints+1; v(i)=0; end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 ttl=['Potential Well(No pertubation), m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 title(ttl);



figure();
for i=1:n
    plot(x,sqrt(2/length)*sin(-i*pi*x/length),s(i)); %wave function
    en(i)=(hbar*i*pi)^2/(2*me*length);
    hold on
end
 tt2=['Potential Well(No pertubation),Analytical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Wave function');
 title(tt2);

en

%numerical solution 
%CASE 2
[energy,phi]=solve_schM(10,npoints,v,1,n);		%call solve_schM
figure()
for i=1:4
       
    
    %j=1+mod(i,7);
      plot(x,(phi(:,i),s1(i));											%plot eigenfunctions
       hold on;
       sprintf(['eigenenergy (',num2str(i),') = ',num2str(energy(i)),' eV'])		%energy eigenvalues
end
 tt3=['Potential Well(No pertubation),Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Wave function');
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
 tt4=['Potential Well(CASE 2), m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 title(tt4);
 
 %numerical solution 

[energy,phi]=solve_schM(10,npoints,v,1,n);		%call solve_schM
figure()
for i=1:4
       
    
    %j=1+mod(i,7);
      plot(x,phi(:,i),s1(i));											%plot eigenfunctions
       hold on;
       sprintf(['eigenenergy CASE 2(',num2str(i),') = ',num2str(energy(i)),' eV'])		%energy eigenvalues
end
 tt4=['Potential Well(CASE 2),Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Wave function');
 title(tt4);
 
 
