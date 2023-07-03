%% Verification of Classic Action using Toffoli Gates

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

X = [0 1;
     1 0];

% Here we are redefining the Toffolli as a function of the control qubits
% c1 c2 and target qubit qt. 

TOFFOLI_ = @(c1c2,t) all(c1c2==k11)*X*t + ~all(c1c2==k11)*t;


%% AND GATE for AB = 00, 01, 10, 11

disp('For AB=00')
disp(TOFFOLI_(k00,ket0))
disp('For AB=01')
disp(TOFFOLI_(k01,ket0))
disp('For AB=10')
disp(TOFFOLI_(k10,ket0))
disp('For AB=11')
disp(TOFFOLI_(k11,ket0))

%Note: The results would be printed in the form of matrix representation of
%|0> and |1> and are hence verified. 

%% NAND GATE for AB = 00, 01, 10, 11

disp('For AB=00')
disp(TOFFOLI_(k00,ket1))
disp('For AB=01')
disp(TOFFOLI_(k01,ket1))
disp('For AB=10')
disp(TOFFOLI_(k10,ket1))
disp('For AB=11')
disp(TOFFOLI_(k11,ket1))

%Note: The results would be printed in the form of matrix representation of
%|0> and |1> and are hence verified. 

%% NOT Gate for A = 0, 1

disp('For A=0')
disp(TOFFOLI_(k11,ket0))
disp('For A=1')
disp(TOFFOLI_(k11,ket1))

%% OR Gate for AB = 00, 01, 10, 11

X2 = kron(X,X); %Since here we have the X Gate acting on both Qubits

disp('For AB=00')
disp(TOFFOLI_(X2*k00,ket1))
disp('For AB=01')
disp(TOFFOLI_(X2*k01,ket1))
disp('For AB=10')
disp(TOFFOLI_(X2*k10,ket1))
disp('For AB=11')
disp(TOFFOLI_(X2*k11,ket1))

%% XOR (CNOT) Gate for AB = 00, 01, 10, 11

disp('For AB=00')
disp(TOFFOLI_(kron(ket1,ket0),ket0))
disp('For AB=01')
disp(TOFFOLI_(kron(ket1,ket0),ket1))
disp('For AB=10')
disp(TOFFOLI_(kron(ket1,ket1),ket0))
disp('For AB=11')
disp(TOFFOLI_(kron(ket1,ket1),ket1))

%%
% Hence, we can conclude that Toffoli Gates can be used to obtain any
% classical operation. We can use Toffoli Gates to develop Half and Full
% Adders for Quantum Circuits. 
