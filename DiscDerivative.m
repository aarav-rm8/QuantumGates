%% AHO using discrete derivative

%Here, I will try to use a well known approximation to obtain the
%approximate Hamiltonian for any 1D Potential. Since the LC-Oscillator and
%Josephson Junction Oscillators can be mapped to 1D systems, we can use
%this method. 


%This method requires 1 assumption however, that the particle is in a
%bounded region. However, if we let N and L be as large as possible, we can
%end up with a fairly close approximation to an unbounded system for lower
%levels of N. Hence in a system of N=1000, we will examine for first 20
%energy levels. 


%% Part 1 : QHO
f1 = figure();
figure(f1)
N = 10000;
dy=1/N;
y = linspace(0,1,N+1);
L = 1000;
m=1;
k=1;

mL2V = m*L^2 * 0.5*k*(y-0.5).^2;



dmain = dy^(-2) + mL2V(2:N);
doff = ones(1,N-2)*-0.5/(dy^2);

mL2H = diag(dmain)+ diag(doff,1) + diag(doff,-1);
H = mL2H / (m*L^2);

E = eig(H);
V;D = eig(H);
subplot(2,2,1)
plot(y,mL2V/(m*L^2))
hold on
yline(E(1:20))
title('Energy Levels (compare to HO)')
xlabel('y')
ylabel('E')

subplot(2,2,2)
bar((0:1:19),E(1:20))
title('Energy Levels (bar graph)')
xlabel('N')
ylabel('Energy')


Lines = [];
for i=1:19
    Lines = [Lines , E(i+1)-E(i)];
end
subplot(2,2,3)
bar(Lines)
title('Transition Energies for Subsequent Levels')
xlabel('N')
ylabel('Energy')

subplot(2,2,4)
xline(0)
yline(0)

l=[]
for num = 1:5

    plot(y(1:N-1),V(:,num))
    l = [l ,'n='+string(num-1)]
    hold on
end
title('First 5 Eigenstates')
xlabel('y')
legend(l)
ylabel('Ψ')
%% Part 2: Recreating the same for a potential with an additional anharmonic term 
f2 = figure();
figure(f2)

k_=1.2;

mL2V = m*L^2 * (0.5*k*(y-0.5).^2 - k_*(y-0.5).^4)

dmain = dy^(-2) + mL2V(2:N);
doff = ones(1,N-2)*-0.5/(dy^2);

mL2H = diag(dmain)+ diag(doff,1) + diag(doff,-1);
H = mL2H / (m*L^2);

E = eig(H);
V;D = eig(H);
subplot(2,2,1)
plot(y,mL2V/(m*L^2))
hold on
yline(E(1:40))
title('Energy Levels (compare to AHO)')
xlabel('y')
ylabel('E')

subplot(2,2,2)
bar((0:1:39),E(1:40))
title('Energy Levels (bar graph)')
xlabel('N')
ylabel('Energy')


Lines = [];
for i=1:39
    Lines = [Lines , E(i+1)-E(i)];
end
subplot(2,2,3)
bar(Lines)
title('Transition Energies for Subsequent Levels')
xlabel('N')
ylabel('Energy')

subplot(2,2,4)
xline(0)
yline(0)

l=[]
for num = 1:5

    plot(y(1:N-1),V(:,num))
    l = [l ,'n='+string(num-1)]
    hold on
end
title('First 5 Eigenstates')
xlabel('y')
legend(l)
ylabel('Ψ')