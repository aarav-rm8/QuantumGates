%%

%% Generation of Toffoli GATE using CNOT, H, T and T† Gates
%Here, we shall attempt to construct the CCX Gate using CNOT,H,T and T†
%Gates. We can use that fact that T† = inverse of T

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

CNOT = [1 0 0 0;
        0 1 0 0;
        0 0 0 1
        0 0 1 0];
I = [1 0; 0 1];
H = (1/sqrt(2))*[1 1;
                1 -1];

% Since it is a 3-Qubit System, we shall define CNOT as a 3-Qubit Gate with
% different sets of Ctrl and Target Qubits used. 

CNOT_12 = kron(CNOT,I);
CNOT_23 = kron(I,CNOT);

CNOT_13 = [1 0 0 0 0 0 0 0;
           0 1 0 0 0 0 0 0
           0 0 1 0 0 0 0 0
           0 0 0 1 0 0 0 0
           0 0 0 0 0 1 0 0
           0 0 0 0 1 0 0 0
           0 0 0 0 0 0 0 1
           0 0 0 0 0 0 1 0]; %This is derivable from the solutions.

T = [1 0;
    0 exp(1i*pi/4)];
Td = inv(T);

%Now we shall multiply the matrices in appropriate order. Will do it
%stepwise and display the matrix after the final operation is done

TOFFOLI_Const = kron(kron(I,I),H); %Step 1
TOFFOLI_Const = CNOT_23*TOFFOLI_Const; %Step 2
TOFFOLI_Const = kron(kron(I,I),Td)*TOFFOLI_Const; %Step 3
TOFFOLI_Const = CNOT_13*TOFFOLI_Const; %Step 4
TOFFOLI_Const = kron(kron(I,I),T)*TOFFOLI_Const; %Step 5
TOFFOLI_Const = CNOT_23*TOFFOLI_Const; %Step 6
TOFFOLI_Const = kron(kron(I,I),Td)*TOFFOLI_Const; %Step 7
TOFFOLI_Const = CNOT_13*TOFFOLI_Const; %Step 8
TOFFOLI_Const = kron(kron(I,T),T)*TOFFOLI_Const; %Step 9
TOFFOLI_Const = kron(CNOT,H)*TOFFOLI_Const; %Step 10
TOFFOLI_Const = kron(kron(T,Td),I)*TOFFOLI_Const; %Step 11
TOFFOLI_Const = CNOT_12*TOFFOLI_Const ;%Step 12

disp('The FINAL MAtrix is as follows:')
disp(round(TOFFOLI_Const,10))
%%
%Hence, it is verified that we obtain the original Toffoli Matrix. 



