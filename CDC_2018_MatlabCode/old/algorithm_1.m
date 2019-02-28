close all; clearvars; clc;

d = 5;

%     mu1         mud
mus = [0.90, 0.95, 1.1, 1.1, 1.2];

Ls = [1 1 2 2 1];
Us = [6 4 4 6 6];

beta = mus(1)^Us(1) * mus(2)^Us(2) * mus(3)^Ls(3) * mus(4)^Ls(4) * mus(5)^Ls(5); % must be < 1, modify for new d!

epsilon = beta + .2; % must be [beta, 1)

ks = Us; j = d; 

while(true)
    if (ks * log(mus)' <= log(epsilon))
        ksStar = ks;
        break;
    else
        if ks(j) == Ls(j), j = j - 1; end
        ks(j) = ks(j) - 1;
    end
end

cvx_begin
    variables kscvx(1,d)
    maximize(sum(kscvx))
    subject to
        kscvx >= Ls;
        kscvx <= Us;
        kscvx * log(mus)' <= log(epsilon);
cvx_end
% floor(sum(kscvx)) should be equal to sum(ks)



