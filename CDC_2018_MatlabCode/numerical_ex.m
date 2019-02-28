%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR: Margaret P. Chapman
% DATE: 12/8/2017
% PURPOSE: Numerical examples for CDC 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DYNAMICS & x0
close all; clearvars; clc;

load('model_id_K14only_aug6.mat'); % load model id, 2 phenotypes
A_Comb = A_STAR{4}; A_BEZ = A_STAR{3}; A_Tram = A_STAR{2};

A_BAR = cell(3,1); A_BAR{1} = A_Comb(1:end-1, 1:end-1); A_BAR{2} = A_BEZ(1:end-1, 1:end-1); A_BAR{3} = A_Tram(1:end-1, 1:end-1);
% DRUG 1 = COMBO, DRUG 2 = BEZ235, DRUG 3 = TRAMETINIB

M = GetCellCtsMATRIX( GetCellCtsCELL( GetK14(P{1}), d{1}, EVEN_DEATH, 1:N_WELL ), 1:N_WELL ); % 15-well data, DMSO, even death, 2 phenotypes

x0 = mean( M(1:end-1, 1:N_WELL), 2); % average [# K14+ live, # K14- live] at time 0

%% GET INDUCED MATRIX 1-NORM (=max eigenvalue for our matrices)

mu1 = norm( A_BAR{1}, 1 ); mu2 = norm( A_BAR{2}, 1 ); mu3 = norm( A_BAR{3}, 1 );

%% SET CONSTANTS

L1 = 2; U1 = 4;             L2 = 2; U2 = 8;             L3 = 2; U3 = 6;

beta = mu1^U1 * mu2^L2 * mu3^L3; % Is this < 1?

if beta < 1, epsilon = (1-beta)/2 + beta; else display('beta not > 1!'); end

epsilon = 0.95; % epsilon \in [beta, 1]

eta = 0.96; % eta > epsilon

%% GET MAXIMAL WAITING TIMES

% drug 1, drug 3, drug 2
ks = [U1 U3 U2]; mus = [mu1 mu3 mu2]; Ls = [L1 L3 L2]; j = 3; % #drugs
                       %small->big
while(true)
    if (ks * log(mus)' <= log(epsilon))
        ksStar = ks;
        break;
    else
        if ks(j) == Ls(j), j = j - 1; end
        ks(j) = ks(j) - 1;
    end
end

%ksStar = [drug 1, drug 3, drug 2]
k1 = ksStar(1); k3 = ksStar(2); k2 = ksStar(3);
  
%% ILLUSTRATE THM1, THM2
                                     % maximal waiting times drug 1, drug 2, drug 3
t00 = 0; x00 = x0; ts = []; ys = []; k = [k1, k2, k3]; ts00 = []; ys00 = []; errs00 = []; errs = [];

% Generate values from the uniform distribution on the interval [a, b]. r = a + (b-a).*rand(100,1); 
a = .9; b = 1.5;

%set to 1 for THM2, set to 0 for THM1
yeserr = 1;

N_Cycle = 40; % # cycles
numloop = 1;
for m = 1 : N_Cycle
    
  p = randperm(3); % P = randperm(N) returns a vector containing a random permutation of the integers 1:N.
  
  while(true)      % try until multiplicative error over entire cycle is bounded by 1/eta
  
    err1x = 1; if yeserr, err1 = a + (b-a).*rand( k(p(1)), 1 )'; for j = 1 : k(p(1)), err1x = err1x * err1(j); end; end;
  
    err2x = 1; if yeserr, err2 = a + (b-a).*rand( k(p(2)), 1 )'; for j = 1 : k(p(2)), err2x = err2x * err2(j); end; end;
  
    err3x = 1; if yeserr, err3 = a + (b-a).*rand( k(p(3)), 1 )'; for j = 1 : k(p(3)), err3x = err3x * err3(j); end; end;
    
    e00 = err1x*err2x*err3x;
    
    if e00 <= 1/eta, break; end 
    
    numloop = numloop+1;
    
  end
    
  xp1 = err1x * A_BAR{p(1)}^k(p(1)) * x00; % apply 1st drug in p
  
  xp1_PLUS_xp2 = err2x * A_BAR{p(2)}^k(p(2)) * xp1; % apply 2st drug in p
  
  xp1_PLUS_xp2_PLUS_xp3 = err3x * A_BAR{p(3)}^k(p(3)) * xp1_PLUS_xp2; % apply 3rd drug in p
  
  ts00 = [ ts00, t00 ]; ys00 = [ ys00, norm( x00, 1 ) ]; % at start of cycle
  
  errs00 = [ errs00, e00 ]; % multiplicative error for current cycle
  
  if yeserr, errs = [ errs, err1, err2, err3 ]; end;
  
  ts = [ ts, t00,            t00 + k(p(1)),  t00 + k(p(1)) + k(p(2)) ];
  
  ys = [ ys, norm( x00, 1 ), norm( xp1, 1 ), norm( xp1_PLUS_xp2, 1 ) ];
  
  x00 = xp1_PLUS_xp2_PLUS_xp3;
  
  t00 =  t00 + k(p(1)) + k(p(2)) + k(p(3));
  
end

ts = [ ts, t00 ]; % store last time point
ys = [ ys, norm( x00, 1 ) ];

%% PLOT RESULTS
myticklabels_days = {'0','','','','','20','','','','','40','','','','','60','','','','','80','','','','','100','','','','','120','','','','','140','','','','','160',};
myticklabels_cycl = {'0','','','','','5','','','','','10','','','','','15','','','','','20','','','','','25','','','','','30','','','','','35','','','','','40',};

figure; FigureSettings;

subplot(3,1,1)
%plot( ts00/2, ys00/norm(x0,1), '*m', 'linewidth', 1 ); hold on
plot( ts/2, ys/norm(x0,1), 'k' );
ylabel('Normalized pop. size'); xlabel('time (days)'); %legend('start of cycle');
xticks([ts00/2, max(ts/2)]); grid on; xticklabels(myticklabels_days);   
                % include tick at 160
if yeserr
subplot(3,1,2)
plot( (0:length(errs)-1)/2, errs, 'ok', 'linewidth', 1 ); %first error starts at time 0
ylabel( 'Error' ); xlabel('time (days)');
xticks([ts00/2, max(ts/2)]); grid on; xticklabels(myticklabels_days);

cycle_time = 1 : N_Cycle; % tick 1 is end of cycle 1,..., tick 40 is end of cycle 40
subplot(3,1,3)
plot( cycle_time, 1/eta*ones(size(cycle_time)), 'xk', 'linewidth', 1); hold on;
plot( cycle_time, errs00, '^b', 'linewidth', 1 );
ylabel('Error product'); xlabel('time (cycles)');
legend('1/\eta'); xticks([0, cycle_time]); grid on; xticklabels(myticklabels_cycl);
        % include tick at 0
end

%% ILLUSRATE SETTLING TIME
% DRUG 1 = COMBO, DRUG 2 = BEZ235, DRUG 3 = TRAMETINIB

gamma = 1/10;

B = b^(k1+k2+k3) * mu2^k2 * mu3^k3; % is B > gamma?

m_gamma = ceil( (( log(gamma) - log(B) )/( log(epsilon) - log(eta) )) + 1 );
    
%% ILLUSTRATE THM1&COR1 OR THM3 - OLD VERSION
% 
% t00 = 0; x00 = x0; ts = []; ys = [];
% 
% %set to 1 for THM3, set to 0 for THM1&COR1
% yeserr = 1;
% 
% for m = 1 : 100
%       
%   err1x = 1; if yeserr, err1 = a + (b-a).*rand( k1, 1 )'; for j = 1 : k1, err1x = err1x * err1(j); end; end;
%   
%   err2x = 1; if yeserr, err2 = a + (b-a).*rand( k2, 1 )'; for j = 1 : k2, err2x = err2x * err2(j); end; end;
%   
%   err3x = 1; if yeserr, err3 = a + (b-a).*rand( k3, 1 )'; for j = 1 : k3, err3x = err3x * err3(j); end; end;
% 
%   xk1 = err1x * A_BAR{1}^k1 * x00;
%   
%   xk1_PLUS_xk2 = err2x * A_BAR{2}^k2 * xk1;
%   
%   xk1_PLUS_xk2_PLUS_xk3 = err3x * A_BAR{3}^k3 * xk1_PLUS_xk2;
%   
%   ts = [ ts, t00,            t00 + k1,       t00 + k1 + k2           ];
%   
%   ys = [ ys, norm( x00, 1 ), norm( xk1, 1 ), norm( xk1_PLUS_xk2, 1 ) ];
%     
%   x00 = xk1_PLUS_xk2_PLUS_xk3;
%   
%   t00 = t00 + k1 + k2 + k3;
%   
% end
% 
% figure; FigureSettings;
% 
% plot( ts/2, ys/norm(x0,1), 'k', 'linewidth', 1 ); ylabel('Normalized # live cancer cells'); xlabel('Time (days)'); 
%     
%     
    
    
