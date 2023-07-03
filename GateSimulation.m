%% Simulation of Single and Multi-Qubit Gates

%In this section, I attempt to simulate the action of gates. I shall be
%using matrix representation of the gates. 

ket0 = [1;0];
ket1 = [0;1];

syms a
syms b
psi = a*ket0 + b*ket1
%I have defined a and b as symbolic variables, making advantage of the
%Symbolic Math Toolbox of MATLAB. The reason for the same is to esure that
%we can truely understand gate action on arbitary state.

format shortG

%% Single Qubit Gate Action

%First, we shall be simulating gate action on single qubits. 

%% 1. Pauli-X
X = [0 1;
    1 0];

disp('Action of X gate on state |0>:')
disp(X*ket0)

disp('Action of X gate on state |1>:')
disp(X*ket1)

disp('Action of X gate on state a|0> + b|1>:')
disp(X*psi)

%% 2. Pauli-Y
Y = [0 -1i;
    1i 0];

disp('Action of Y gate on state |0>:')
disp(Y*ket0)

disp('Action of Y gate on state |1>:')
disp(Y*ket1)

disp('Action of Y gate on state a|0> + b|1>:')
disp(Y*psi)

%% 3. Pauli-Z
Z = [1 0;
    0 -1];

disp('Action of Z gate on state |0>:')
disp(Z*ket0)

disp('Action of Z gate on state |1>:')
disp(Z*ket1)

disp('Action of Z gate on state a|0> + b|1>:')
disp(Z*psi)

%% 4. Hadamard Gate
H = (1/sqrt(2))*[1 1;
                1 -1];

disp('Action of H gate on state |0>:')
disp(H*ket0)

disp('Action of H gate on state |1>:')
disp(H*ket1)

disp('Action of H gate on state a|0> + b|1>:')
disp(H*psi)

%% 5. S Gate
S = [1 0;
    0 j];

disp('Action of S gate on state |0>:')
disp(S*ket0)

disp('Action of S gate on state |1>:')
disp(S*ket1)

disp('Action of S gate on state a|0> + b|1>:')
disp(S*psi)


%% 6. T Gate
T = [1 0;
    0 exp(i*pi/4)];

disp('Action of T gate on state |0>:')
disp(T*ket0)

disp('Action of T gate on state |1>:')
disp(T*ket1)

disp('Action of T gate on state a|0> + b|1>:')
disp(T*psi)

%% Verification of basic gate identites

HZH = H*Z*H
HXH = H*X*H
disp('(X+Z)/√2=')
disp((X+Z)/sqrt(2))
XZ = X*Z
ZX = Z*X
iY = i*Y

%% Action of 2-Qubit Gates

%For 2-qubit gates, we have four basis states which stem from tensor 
% product of the basis states of a single qubit

syms c
syms d

k00 = kron(ket0,ket0);
k01 = kron(ket0,ket1);
k10 = kron(ket1,ket0);
k11 = kron(ket1,ket1);

psi_ = a*k00+b*k01+c*k10+d*k11


%% CNOT-Gate

CNOT = [1 0 0 0;
        0 1 0 0;
        0 0 0 1
        0 0 1 0];

disp('Action of CNOT gate on state |00>:')
disp(CNOT*k00)

disp('Action of CNOT gate on state |01>:')
disp(CNOT*k01)

disp('Action of CNOT gate on state |10>:')
disp(CNOT*k10)

disp('Action of CNOT gate on state |11>:')
disp(CNOT*k11)

disp('Action of CNOT gate on state a|00> + b|01> + c|10> + d|11>:')
disp(CNOT*psi_)


%% CZ Gate
CZ = [1 0 0 0;
      0 1 0 0;
      0 0 1 0
      0 0 0 -1];

disp('Action of CZ gate on state |00>:')
disp(CZ*k00)

disp('Action of CZ gate on state |01>:')
disp(CZ*k01)

disp('Action of CZ gate on state |10>:')
disp(CZ*k10)

disp('Action of CZ gate on state |11>:')
disp(CZ*k11)

disp('Action of CZ gate on state a|00> + b|01> + c|10> + d|11>:')
disp(CZ*psi_)


%% SWAP GATE
SWAP = [1 0 0 0;
        0 0 1 0;
        0 1 0 0
        0 0 0 1];

disp('Action of SWAP gate on state |00>:')
disp(SWAP*k00)

disp('Action of SWAP gate on state |01>:')
disp(SWAP*k01)

disp('Action of SWAP gate on state |10>:')
disp(SWAP*k10)

disp('Action of SWAP gate on state |11>:')
disp(SWAP*k11)

disp('Action of SWAP gate on state a|00> + b|01> + c|10> + d|11>:')
disp(SWAP*psi_)

%% 3 Qubit Gate: Toffoli Gate

syms e
syms f
syms g
syms h
%first, we shall construct our 3-qubit basis states
k000 = kron(ket0,k00);
k001 = kron(ket0,k01);
k010 = kron(ket0,k10);
k011 = kron(ket0,k11);
k100 = kron(ket1,k00);
k101 = kron(ket1,k01);
k110 = kron(ket1,k10);
k111 = kron(ket1,k11);

psi__ = a*k000 + b*k001 + c*k010 + d*k011 + e*k100 + f*k101 + g*k110 + h*k111

TOFFOLI = [1 0 0 0 0 0 0 0;
           0 1 0 0 0 0 0 0
           0 0 1 0 0 0 0 0
           0 0 0 1 0 0 0 0
           0 0 0 0 1 0 0 0
           0 0 0 0 0 1 0 0
           0 0 0 0 0 0 0 1
           0 0 0 0 0 0 1 0];

% We shall not be testing out the Toffoli for individual cases, but we
% shall instead be checking it for our general state only and make our
% inferences from there. It shall be verified if g and h would swap in
% the final statevector.

disp('Action of TOFFOLI gate on general state a|000> + b|001> + c|010> + d|011>')
disp('+ e|100> + f|101> + g|110> + h|111>')
disp(TOFFOLI*psi__)


