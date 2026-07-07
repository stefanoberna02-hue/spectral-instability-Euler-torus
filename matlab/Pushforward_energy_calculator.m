%% tg_fourier_generic.m
clear; clc;

m = 2;
n = 2;

% Each row [a,j,k] represents a*exp(1i*(j*x+k*y)).
% We choose g0 = cos(3x+y), written in complex Fourier modes.
g0 = [
    1/2,  3,  1;
    1/2, -3, -1
];

% Since J is skew-adjoint, J^* = -J.
% Hence f = -Jg0 belongs to Range(J^*) = Ker(J)^\perp.
f = -applyJmn(g0,m,n);

Hf = applyHmn(f,m,n);

energy = innerL2(f,Hf);

disp('Initial function g0 = cos(3x+y):');
disp(g0);

disp('Function f = -J g0:');
disp(f);

disp('H f:');
disp(Hf);

disp('<f,Hf>_L2 =');
disp(energy);

if real(energy) < 0
    disp('Negative energy found in Ker(J)^\perp.');
else
    disp('This initial choice does not give negative energy.');
end


%% Functions

function HF = applyHmn(F,m,n)
% Applies H_{m,n} = I + (m^2+n^2) Delta^{-1}.

    HF = zeros(size(F));

    for r = 1:size(F,1)
        a = F(r,1);
        j = F(r,2);
        k = F(r,3);

        R2 = j^2 + k^2;

        if R2 == 0
            error('H_{m,n} is not defined on the zero Fourier mode.');
        end

        multiplier = 1 - (m^2+n^2)/R2;

        HF(r,:) = [a*multiplier, j, k];
    end
end


function JF = applyJmn(F,m,n)
% Applies J_{m,n} = -{psi_{m,n},.}, psi_{m,n}=-sin(mx)sin(ny).

    JF = zeros(4*size(F,1),3);
    row = 1;

    for r = 1:size(F,1)
        a = F(r,1);
        j = F(r,2);
        k = F(r,3);

        JF(row,:) = [a*(m*k - n*j)/4,  j+m, k+n];
        row = row + 1;

        JF(row,:) = [a*(-m*k - n*j)/4, j+m, k-n];
        row = row + 1;

        JF(row,:) = [a*(m*k + n*j)/4,  j-m, k+n];
        row = row + 1;

        JF(row,:) = [a*(n*j - m*k)/4,  j-m, k-n];
        row = row + 1;
    end
end


function val = innerL2(F,G)
% Standard complex L2 product on [0,2*pi]^2:
% <f,g> = integral conj(f) g.

    val = 0;

    for r = 1:size(F,1)
        a = F(r,1);
        j = F(r,2);
        k = F(r,3);

        for s = 1:size(G,1)
            b = G(s,1);
            p = G(s,2);
            q = G(s,3);

            if j == p && k == q
                val = val + conj(a)*b*4*pi^2;
            end
        end
    end
end