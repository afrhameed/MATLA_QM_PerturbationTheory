npoints=200;
me=9.10938188e-31; %mass of electron 
length=10e-9;   %10nm
x=0:length/npoints:length;

hbar=1.054571596e-34;	%Planck's constant (x10^34 J s)
k=(0.0038/0.0602)*1.0e28;
alpha=(0.001570731731182)^2/(4.934396342684429e+04) ;

n=4  % up to first 4 bound states
s1=char('.y','.k','.r','.g','.b','.m','.c');
s=char('b','r','y','m','b','m','c');								%plot curves in different colors

%CASE 1___________________________________

figure();
for i=1:n
    plot(x,((sqrt(2/length)*sin(-i*pi*x/length)).^2)*alpha,s(i)); %wave function
    en(i)=(hbar*i*pi)^2/(2*me*length);
    hold on
end
 tt2=['Case 1,Analytical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Probability Density');
 title(tt2);

en.*k

%CASE 2___________________________________

figure();
for i=1:n
    plot(x,((sqrt(2/length)*sin(-i*pi*x/length)).^2)*alpha,s(i)); %wave function
    en(i)=(hbar*i*pi)^2/(2*me*length);
    hold on
end
 tt2=['Case 1,Analytical sol. m* = ',num2str(me),'m0, Length = ',num2str(length),'nm'];
 legend('n=1','n=2','n=3','n=4');
 xlabel('Distance (nm)'),ylabel('Probability Density');
 title(tt2);

en.*k

