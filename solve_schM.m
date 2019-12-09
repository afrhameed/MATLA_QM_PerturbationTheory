function [e,phi]=solve_schM(length,n,v,mass,num_sol)
% Solve Schrodinger equation for energy eigenvalues and eigenfunctions.
% [energy, phi_fun]=solve_sch(length,number,potential,mass,sol_num)
% solves time-independent Schrodinger equation in region bounded by 0<=x<=length.
% Potential is infinite outside this region.
%
% length			length of region (nm)
% n				number of sample points 
% v				potential inside region (eV)
% mass			effective electron mass
% sol_num		number of solutions sought
%
% e				energy eigenvalue (eV)
% phi				eigenfunction with eigenvalue = e
%

hbar=1.054571596;										%Planck's constant (x10^34 J s)
hbar2=hbar^2;
echarge=1.602176462;									%electron charge (x10^19 C)
baremass=9.10938188;									%bare electron mass (x10^31 kg)

const=hbar2/baremass/echarge;						%(hbar^2)/(echarge*1nm^2*m0)
const=const/mass;										%/m-effective
deltax=length/n;										%x-increment = length(nm)/n
deltax2=deltax^2;
const=const/deltax2;

for i=2:n
   d(i-1)=v(i)+const;								%diagonal matrix element				
   offd(i-1)=const/2;								%off-diagonal matrix element
end

t1=-offd(2:n-1);							
Hmatrix=diag(t1,1);									%upper diagonal of Hamiltonian matrix							
Hmatrix=Hmatrix+diag(t1,-1);						%add lower diagonal of Hamiltonian matrix

t2=d(1:n-1);
Hmatrix=Hmatrix+diag(t2,0);						%add diagonal of Hamiltonian matrix

[phi,te]=eigs(Hmatrix,num_sol,'SM');			%use matlab function eigs to find num_sol eigenfunctions and eigenvalues

for i=1:num_sol
   e(i)=te(i,i);										%return energy eigenvalues in vector e
end
phi=[zeros(1,num_sol);phi;zeros(1,num_sol)]; %wave function is zero at x = 0 and x = length
return

