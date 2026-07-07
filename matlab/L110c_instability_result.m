%% L110c_instability_result.m
clear; clc;

% Assumption: each row [a,j,k] represents a*cos(j*x+k*y).
% We work on the standard torus [0,2*pi]^2 with the standard L2 product.

m = 2;
n = 2;

% Successful function g leading to the instability result.
% The function belongs to the odd-odd cosine class relevant to the (2,2)
% case. Applying J_{2,2}^* = -J_{2,2} gives
%
%     f = J_{2,2}^* g = -J_{2,2} g.
%
% Hence f belongs to Range(J_{2,2}^*), which is contained in
% Ker(J_{2,2})^\perp. The resulting vector has negative H_{2,2}-energy.
g = [
   -1,  3, -1;
    1, -1,  3
];

% Since J is skew-adjoint, J^* = -J.
f = scaleModes(applyJ(g,m,n), -1);

Hf = applyH(f,m,n);

energy = innerL2(f,Hf);

disp('f = -J g is represented by:');
disp(f);

disp('<f, H f> =');
disp(energy);

if energy < 0
    disp('Negative energy found in Ker(J)^\perp.');
else
    disp('The chosen test function does not give negative energy.');
end


%% Local functions

function out = applyH(F,m,n)
% Applies H_{m,n} = Id + (m^2+n^2) Delta^{-1}.

    out = F;

    for r = 1:size(F,1)
        j = F(r,2);
        k = F(r,3);
        R2 = j^2 + k^2;

        if R2 == 0
            error('H is not defined on the zero Fourier mode.');
        end

        multiplier = 1 - (m^2+n^2)/R2;
        out(r,1) = F(r,1) * multiplier;
    end

    out = collectModes(out);
end


function out = applyJ(F,m,n)
% Applies J_{m,n}f = -{psi_{m,n},f}, with
% psi_{m,n} = -sin(mx)sin(ny).

    out = zeros(4*size(F,1),3);
    row = 1;

    for r = 1:size(F,1)
        a = F(r,1);
        j = F(r,2);
        k = F(r,3);

        out(row:row+3,:) = [
            a*(m*k - n*j)/4,   j+m, k+n;
           -a*(m*k + n*j)/4,   j+m, k-n;
            a*(m*k + n*j)/4,   j-m, k+n;
            a*(n*j - m*k)/4,   j-m, k-n
        ];

        row = row + 4;
    end

    out = collectModes(out);
end


function val = innerL2(F,G)
% Standard L2 inner product on [0,2*pi]^2.
% Each row represents a*cos(jx+ky).

    F = collectModes(F);
    G = collectModes(G);

    val = 0;

    for r = 1:size(F,1)
        a = F(r,1);
        j = F(r,2);
        k = F(r,3);

        for s = 1:size(G,1)
            b = G(s,1);
            p = G(s,2);
            q = G(s,3);

            sameMode = (j == p && k == q);

            if sameMode
                if j == 0 && k == 0
                    normFactor = 4*pi^2;
                else
                    normFactor = 2*pi^2;
                end

                val = val + a*b*normFactor;
            end
        end
    end
end


function out = scaleModes(F,c)
% Multiplies all Fourier coefficients by c.

    out = F;
    out(:,1) = c*out(:,1);
    out = collectModes(out);
end


function out = collectModes(F)
% Combines equal cosine modes, using cos(-jx-ky)=cos(jx+ky).

    tol = 1e-12;

    if isempty(F)
        out = F;
        return;
    end

    % Remove numerically zero coefficients.
    F(abs(F(:,1)) < tol,:) = [];

    if isempty(F)
        out = zeros(0,3);
        return;
    end

    % Canonicalize the cosine symmetry: (j,k) ~ (-j,-k).
    for r = 1:size(F,1)
        j = F(r,2);
        k = F(r,3);

        if j < 0 || (j == 0 && k < 0)
            F(r,2) = -j;
            F(r,3) = -k;
        end
    end

    modes = F(:,2:3);
    coeffs = F(:,1);

    [uniqueModes,~,idx] = unique(modes,'rows');
    newCoeffs = accumarray(idx,coeffs);

    out = [newCoeffs, uniqueModes];

    % Remove modes cancelled by summation.
    out(abs(out(:,1)) < tol,:) = [];
end
