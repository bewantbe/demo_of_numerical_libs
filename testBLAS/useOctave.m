% use Octave (usually link to libatlas.so)
toct = @() fprintf('  time: %5.3f sec.\n', toc());

m = 1000;
n = 1000;

disp('Allocating memory:');
tic();
A = zeros(m,n);
B = zeros(m,n);
C = zeros(m,n);
D = zeros(m,n);
toct();

disp('Memory cost:');
fprintf('  A(%d*%d): %.2fMB\n', m, n, whos('A').bytes/2^20);
fprintf('  B(%d*%d): %.2fMB\n', m, n, whos('B').bytes/2^20);
fprintf('  C(%d*%d): %.2fMB\n', m, n, whos('C').bytes/2^20);
fprintf('  D(%d*%d): %.2fMB\n', m, n, whos('D').bytes/2^20);
fprintf(' total: %.2fMB\n\n', (whos('A').bytes + whos('B').bytes + whos('C').bytes + whos('D').bytes)/2^20);

fprintf('Problem size: %d * %d\n', m, n);
disp('generating: A[i][j]');
tic();
A = floor(rand(m, n) * 2^9);
toct();

disp('copying: B = A');
tic();
B = A;                   % Octave has lazy copy
B(1) = B(1) + 1e-300;    % so here try to cause the copy happen
toct();

disp('calculating: C = A * B');
tic();
C = A*B;
toct();

n_rhs = 1;
fprintf('solving: D: A * D = C  (NRHS = %d)\n', n_rhs);
tic();
D(:,1) = A \ C(:,1);
toct();
disp('cal: error (B-D)_Abs');
fprintf('  mean abs error = %e\n', mean(abs((B-D)(1:n_rhs))));
fprintf('  max  abs error = %e\n', max (abs((B-D)(1:n_rhs))));
fprintf('  abs norm  = %e\n', mean(abs(B(1:n_rhs))));

tic();
n_rhs = m;
fprintf('solving: D: A * D = C  (NRHS = %d)\n', n_rhs);
tic();
D = A \ C;
toct();
disp('cal: error (B-D)_Abs');
fprintf('  mean abs error = %e\n', mean(abs((B-D)(1:n_rhs))));
fprintf('  max  abs error = %e\n', max (abs((B-D)(1:n_rhs))));
fprintf('  abs norm  = %e\n', mean(abs(B(1:n_rhs))));

