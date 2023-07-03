%% Verification of gate compositions via Rotation and Phase Gates 

%% Redefining Basis States and Basic Gates 
format short

ket0 = [1;0];
ket1 = [0;1];
k00 = kron(ket0,ket0);
k01 = kron(ket0,ket1);
k10 = kron(ket1,ket0);
k11 = kron(ket1,ket1);
k000 = kron(ket0,k00);
k001 = kron(ket0,k01);
k010 = kron(ket0,k10);
k011 = kron(ket0,k11);
k100 = kron(ket1,k00);
k101 = kron(ket1,k01);
k110 = kron(ket1,k10);
k111 = kron(ket1,k11);


syms a
syms b
syms c
syms d
syms e
syms f
syms g
syms h
psi = a*ket0 + b*ket1;
psi_ = a*k00+b*k01+c*k10+d*k11;
psi__ = a*k000 + b*k001 + c*k010 + d*k011 + e*k100 + f*k101 + g*k110 + h*k111;

I = diag(ones(1,2));
X = [0 1;
    1 0];
Y = [0 -1i;
    1i 0];
Z = [1 0;
    0 -1];
H = (1/sqrt(2))*[1 1;
                1 -1];
S = [1 0;
    0 1i];
T = [1 0;
    0 exp(1i*pi/4)];
CNOT = [1 0 0 0;
        0 1 0 0;
        0 0 0 1
        0 0 1 0];
CZ = [1 0 0 0;
      0 1 0 0;
      0 0 1 0
      0 0 0 -1];
SWAP = [1 0 0 0;
        0 0 1 0;
        0 1 0 0
        0 0 0 1];
TOFFOLI = [1 0 0 0 0 0 0 0;
           0 1 0 0 0 0 0 0
           0 0 1 0 0 0 0 0
           0 0 0 1 0 0 0 0
           0 0 0 0 1 0 0 0
           0 0 0 0 0 1 0 0
           0 0 0 0 0 0 0 1
           0 0 0 0 0 0 1 0];

%% Defining some additional gates: Rotation and Phase Gates

Rx = @(th) [cos(th/2) , -1i*sin(th/2);
            -1i*sin(th/2), cos(th/2)];
Ry = @(th) [cos(th/2) , -sin(th/2);
            sin(th/2), cos(th/2)];
Rz = @(th) [exp(-1i*th/2) , 0;
            0, exp(+1i*th/2)];

P = @(phi) [1, 0;
           0, exp(i*phi)];

%% We shall obtain our basic gates using Rx(θ), Ry(θ), Rz(θ), P(φ), CNOT
% For this family we will be attempting to obtain the gates through
% rotations and shall be verifying their matrices.
% For our purposes, we shall ignore any additional Global Phase which might
% get added on since it does not have significance for our computations.

% Pauli-X
disp('X = Rx(π) = ')
disp(Rx(pi))
disp('Here, an additional global phase of -1i is obtained')

% Pauli-Y
disp('Y = Ry(π) =')
disp(Ry(pi))
disp('Here, an additional global phase of -1i is obtained')

% Pauli-Z
disp('Z = P(π) =')
disp(P(pi))


%Hadamard
disp('H = Rx(π)Ry(π/2) =')
disp(Rx(pi)*Ry(pi/2))
disp('Here, an additional global phase of -1i is obtained')

%S and T Gates:
disp('S = P(pi/2) = ')
disp(P(pi/2))
disp('T = P(pi/4) = ')
disp(P(pi/4))

% Building CZ using H(=Rx(π)Ry(π/2)) and CNOT

H_ = Rx(pi)*Ry(pi/2);
disp('CZ = I⊗H * CNOT * I⊗H')
CZ_ = kron(I,H_)*CNOT*kron(I,H_)

% Building SWAP using 2 CNOT Gates
%NOTE: The CNOT Gate Matrix for when the ctrl and the target qubit are
%interechanged would be different. The construction of the inverted CNOT
%matrix would be: 
CNOT_Flip = [1 0 0 0;
         0 0 0 1
         0 0 1 0
         0 1 0 0]; 

disp('SWAP = CNOT*CNOT_Flip*CNOT =')
disp(CNOT*CNOT_Flip*CNOT)








