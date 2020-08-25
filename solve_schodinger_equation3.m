clear,clc;
% PARAMETERS:
L = 1; % Interval Length
N = 1000; % No of points
table = readtable('Data_1D_Configurations_HardCore.N145/Crystal_POSCAR_0.99.txt','ReadVariableNames',0);
coords = table.Var1;
coords = rescale(coords);
x = linspace(0,L,N)'; % Coordinate vector
dx = x(2) - x(1); % Coordinate step
% POTENTIAL
w = L/1000;
U = 1;
for i= 1:145
    U = U+ heaviside(x-(coords(i)-w)) - heaviside(x-(coords(i)+w));
end
U = 10000 * U;
% Three-point finite-difference representation of Laplacian
% using sparse matrices, where you save memory by only
% storing non-zero matrix elements
% 3 point Laplacian
 Lap = (-2*diag(ones(N,1),0) + diag(ones((N-1),1),1)+diag(ones((N-1),1),-1))/(dx^2);
 Lap(1,1) = 0; Lap(1,2) = 0; Lap(2,1) = 0; % So that f(0) = 0
Lap(N,N-1) = 0; Lap(N-1,N) = 0; Lap(N,N) = 0;% So that f(L) = 0

% 5 point Laplacian
%Lap = (-30*diag(ones(N, 1), 0) + 16*diag(ones((N-1), 1), 1) + 16*diag(ones((N-1), 1), -1) - diag(ones((N-2), 1), 2) - diag(ones((N-2), 1), -2)) / (12*dx^2);
%29Lap(1, 1:3) = 0; Lap(1:3, 1) = 0; Lap(N, N-2:N) = 0; Lap(N-2:N, N) = 3;

% Total Hamiltonian
hbar = 1; m = 1; % constants for Hamiltonian
H = -1/2*(hbar^2/m)*Lap +spdiags(U,0,N,N);

% Find lowest nmodes eigenvectors and eigenvalues of sparse matrix
nmodes = 3; options.disp = 0;
[V,E] = eig(H); % find eigs
[E,ind] = sort(diag(E));% convert E to vector and sort low to high
V = V(:,ind); % rearrange corresponding eigenvectors
% Generate plot of lowest energy eigenvectors V(x) and U(x)
Usc = U*max(abs(V(:)))/max(abs(U)); % rescale U for plotting
probdensity = abs(V).^2;
plot(x,probdensity(:,3),'r','LineWidth',4);
axis square;
box on;
set(gca,'linewidth',2,'fontsize',20,'fontname','Times');% plot V(x) and rescaled U(x) ,x,Usc,'--k'

%figure;
%V = abs(V).^2;
%Vmean = mean(V);
%V2 = V - Vmean;
%phi = fft(V2);
%yshift = fftshift(phi(:,3));
%f = (0:N-1)/N * (1/0.001) - 500;
%plot(f,abs(yshift));
%plot(x,V(:,3));
