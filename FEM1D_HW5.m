clc; close all; clear all;
% ## START INPUT ##
% #FEM input
% Number of Elements
Nelem = 32;
% Poly order
Poly = 3;

% #Problem input
A = pi*1e-4;
E1=200e9;
E2=70e9;
rho1=7800;
rho2=2700;
L = 1;
P = 320*9.81;
% ## END INPUT ##
if (Poly == 1)
 Nnod = Nelem+1;
elseif (Poly == 2)
 Nnod = 2*Nelem+1;
elseif(Poly == 3)
 Nnod = 3*(Nelem-1)+4;
end
coord = linspace(0,L,Nnod);
EA = zeros(Nelem);
tx = zeros(Nelem);
for i=1:Nelem
 % Get the nodes and EFT
 if(Poly == 1)
 x1=coord(i);
 x2=coord(i+1);
 xNodes(i,:) = [x1 x2];
 EFT(i,:) = [i i+1];
 elseif(Poly == 2)
 x1=coord(i*2-1);
 x2=coord(i*2);
 x3=coord(i*2+1);
 xNodes(i,:) = [x1 x2 x3];
 EFT(i,:) = [2*i-1 2*i 2*i+1];
 elseif(Poly == 3)
 x1=coord((i-1)*3+1);

 x2=coord((i-1)*3+2);
 x3=coord((i-1)*3+3);
 x4=coord((i-1)*3+4);
 xNodes(i,:) = [x1 x2 x3 x4];
 EFT(i,:) = [(i-1)*3+1 (i-1)*3+2 (i-1)*3+3 (i-1)*3+4];
 end
 xi = x1;
 xj = xNodes(i,end);
 % Get material properties and body force (self-weight) that is constant
 % over each element.
 if (xj <= 0.25)
 EA(i) = E1*A;
 tx(i)=7800*9.81*A;
 elseif (xi >= 0.25 && xj <= 0.5)
 EA(i) = E2*A;
 tx(i)=2700*9.81*A;
 elseif(xi >= 0.5 && xj <= 0.75)
 EA(i) = E1*A;
 tx(i)=7800*9.81*A;
 elseif(xi >= 0.75 && xj <= 1.0)
 EA(i) = E2*A;
 tx(i)=2700*9.81*A;
 else
 error('Element cannot lie between different materials');
 end
end
% Assemble the global stiffness matrix K (same as in HW #1)
K = zeros(Nnod);
F = zeros(Nnod,1);
for i=1:Nelem
 % Call the subroutine to compute the element stiffness matrix Ke
 if(Poly==1)
 [Ke,Fe] = ElemK1Dp1(xNodes(i,:),EA(i),tx(i));
 elseif(Poly==2)
 [Ke,Fe] = ElemK1Dp2(xNodes(i,:),EA(i),tx(i));
 elseif(Poly==3)
 [Ke,Fe] = ElemK1Dp3(xNodes(i,:),EA(i),tx(i));
 end
 % Element freedom table for the element
 EFTe = EFT(i,:);
 % Add the element stiffness matrix Ke to the global stiffness matrix
 % using the element freedom table EFTe
 K(EFTe, EFTe) = K(EFTe, EFTe) + Ke;
 F(EFTe,1) = F(EFTe,1) + Fe;
end
%Applying displacement boundary conditions
Khat = K; Fhat = F;
Khat(1, :)=0;
Khat(:, 1)=0;
Khat(1, 1) = 1.0;
Fhat(1) = 0;
%Applying boundary conditions at the end
Fhat(Nnod) = Fhat(Nnod)+P;
% Solve the linear system
U = Khat\Fhat;
% Compute the reaction forces F using K and U
F = K*U;
% Compute the axial force for each element and store in P
Paxial = zeros(Nelem,1);

% Store location to be ploted where epsilon and stress are calculated
pgaussPlot = [];
epsilonPlot = [];
sigmaPlot =[];
for i=1:Nelem
% Element freedom table for the element
EFTe = EFT(i,:);
% Extract the external displacement vector Ue from the global displacement vector
Ue = U(EFTe);
% Compute the epsilon, stress and axial force
if (Poly==1)
 [Pe,epsilon,sigma,gaussInGlobalCoord] = CompIntForce1Dp1(xNodes(i,:), EA(i), Ue,
A);
elseif (Poly==2)
 [Pe,epsilon,sigma,gaussInGlobalCoord] = CompIntForce1Dp2(xNodes(i,:), EA(i), Ue,
A);
elseif (Poly==3)
 [Pe,epsilon,sigma,gaussInGlobalCoord] = CompIntForce1Dp3(xNodes(i,:), EA(i), Ue,
A);
end
pgaussPlot = [pgaussPlot gaussInGlobalCoord];
epsilonPlot = [epsilonPlot epsilon];
sigmaPlot = [sigmaPlot sigma];
end
%get analytical solution
[xx,uu,ee,ss] = PlotExact();
% plot interpolated FEM solution
%interpolation throught the element (tbd - problem with more than 1 element)
if Nelem==1
xxi = linspace(-1,1,1000);
for i=1:1000
 if (Poly==1)
 Nxi = @(xi)[0.5*(1 - xi) , 0.5*(1 + xi)];
 uuFEM(i) = Nxi(xxi(i))*U;
 elseif (Poly==2)
 Nxi = @(xi)[xi * (xi - 1.0) * 0.5 , 1.0 - xi*xi , xi * (xi + 1.0) * 0.5];
 uuFEM(i) = Nxi(xxi(i))*U;
 elseif (Poly==3)
 Nxi = @(xi)[9.0/16. * (1.0/9.0 - xi*xi) * (xi - 1.0) , 27.0/16.0 * (1.0 -
xi*xi) * (1.0/3.0 - xi) , 27.0/16.0 * (1.0 - xi*xi) * (1.0/3.0 + xi) , -9.0/16. *
(1.0/9.0 - xi*xi) * (xi + 1.0)];
 uuFEM(i) = Nxi(xxi(i))*U;
 end
end
xxFEM = linspace(0,L,1000);
end
if Nelem>1
 xxFEM = coord;
 uuFEM = U;
end
%plot displacement
subplot(3,1,1);
plot(xxFEM,uuFEM);
hold on
plot(xx,uu,'--');
title('Displacement');
xlabel('Coodinate');

ylabel('Displacement');
legend('FEM','analytical')
%plot strain
subplot(3,1,2);
plot(pgaussPlot,epsilonPlot,'*')
hold on
plot(xx,ee,'--');
title('Strain');
xlabel('Coodinate');
ylabel('Strain');
legend('FEM','analytical')
%plot stress
subplot(3,1,3);
plot(pgaussPlot,sigmaPlot,'*')
hold on
plot(xx,ss,'--');
title('Stress');
xlabel('Coodinate');
ylabel('Stress');
legend('FEM','analytical')
format long
% Strain Energy
disp('FEM Strain Energy')
0.5*U'*K*U


