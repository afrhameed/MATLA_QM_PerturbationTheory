clear
clf;


npoints=200;
me=9.10938188e-31; %mass of electron 
length=10e-9;   %10nm
x=0:length/npoints:length;

hbar=1.054571596;										%Planck's constant (x10^34 J s)

n=4  % up to first 4 bound states
s1=char('.y','.k','.r','.g','.b','.m','.c');
s=char('b','r','y','m','b','m','c');								%plot curves in different colors

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
 tt4=['Potential Well(CASE 3), m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
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
 tt4=['Potential Well(CASE 3),Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Wave function');
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
 tt4=['Potential Well(CASE 4), m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
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
 tt4=['Potential Well(CASE 4),Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Wave function');
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
 tt4=['Potential Well(CASE 5), m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 title(tt4);
 
 %numerical solution 

[energy,phi]=solve_schM(10,npoints,v,1,n);		%call solve_schM
figure()
for i=1:4
       
    
    %j=1+mod(i,7);
      plot(x,phi(:,i),s1(i));											%plot eigenfunctions
       hold on;
       sprintf(['eigenenergy CASE 5(',num2str(i),') = ',num2str(energy(i)),' eV'])		%energy eigenvalues
end
 tt4=['Potential Well(CASE 2),Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Wave function');
 title(tt4);
 
 %CASE 6____________________________________________________

%potential well diagram
figure();
for i=1:npoints+1
        v(i)=0.01*sin(i);
    end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 tt4=['Potential Well(CASE 6), m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
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
 tt4=['Potential Well(CASE 6),Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Wave function');
 title(tt4);
 
 
  %CASE 7____________________________________________________

%potential well diagram
figure();
for i=1:npoints+1
        v(i)=0.01*sin(i/20);
    end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 tt4=['Potential Well(CASE 7), m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
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
 tt4=['Potential Well(CASE 7),Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Wave function');
 title(tt4);
 
  %CASE 8____________________________________________________

%potential well diagram
figure();
for i=1:npoints+1
        v(i)=0.01*cos(i/20);
    end						    %potential (eV)
 plot(x,v,'b');xlabel('Distance (nm)'),ylabel('Potential energy, (eV)');
 tt4=['Potential Well(CASE 8), m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
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
 tt4=['Potential Well(CASE 8),Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Wave function');
 title(tt4);
 
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
 tt4=['Potential Well(CASE 9), m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
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
 tt4=['Potential Well(CASE 9),Numerical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Wave function');
 title(tt4);