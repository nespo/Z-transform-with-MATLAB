%% Week 3 – Z-Transform (A1–A5)
% Author: Md Al-Naim
% How to use:
%   - Put this file in a folder, open MATLAB there, and run it (F5).
%   - All images/text outputs will be saved to ./figures
%
% Toolboxes:
%   - Signal Processing Toolbox is used for zplane/freqz.

clear; close all; clc;

% ------------------------ IO setup ------------------------
outdir = fullfile(pwd,'figures');
if ~exist(outdir,'dir'), mkdir(outdir); end

% Make DSP unit-step convention explicit: u[0] = 1
sympref('HeavisideAtOrigin',1);

%% ===================== A1: Finite sequences ============================
% (a) i) x[n] = {1,2,5} at n={0,1,2}  ->  X(z) = 1 + 2 z^{-1} + 5 z^{-2}
syms z n
X1 = 1 + 2*z^(-1) + 5*z^(-2);
disp('A1(i)  X1(z) ='); pretty(X1)

%     ii) x[n] = {0,3,0,4} at n={0,1,2,3} -> X(z) = 3 z^{-1} + 4 z^{-3}
X2 = 0 + 3*z^(-1) + 0*z^(-2) + 4*z^(-3);
disp('A1(ii) X2(z) ='); pretty(X2)

% Save a tiny text summary
fid = fopen(fullfile(outdir,'A1_summary.txt'),'w');
fprintf(fid,'A1(i): X1(z) = 1 + 2 z^{-1} + 5 z^{-2}\n');
fprintf(fid,'A1(ii): X2(z) = 3 z^{-1} + 4 z^{-3}\n');
fprintf(fid,'ROC: finite-length -> entire z-plane except possibly z=0 or z=∞ (no poles)\n');
fclose(fid);

%% ================= A2: Infinite sequences & ROC ========================
% (a) x[n] = a^n u[n], a = 0.6 -> X(z) = 1/(1 - a z^{-1}), ROC |z| > |a|
clear n z; syms n z
assume(n,'integer'); assumeAlso(n>=0);

a = sym(3)/5;                  % exact 0.6
Xa = ztrans(a^n * heaviside(n), n, z);
disp('A2(a) X(z) for a=0.6:'); pretty(Xa)
roc_a = '|z| > 0.6';

% (b) x[n] = (-0.8)^n u[n] -> X(z) = 1/(1 - (-0.8) z^{-1}) = 1/(1 + 0.8 z^{-1})
a2 = -sym(4)/5;                % exact -0.8
Xb = ztrans(a2^n * heaviside(n), n, z);
disp('A2(b) X(z) for a=-0.8:'); pretty(Xb)
roc_b = '|z| > 0.8';

% (c) x[n] = -(0.9)^n u[-n-1]  (left-sided) -> X(z) = -(z)/(z - 0.9), ROC |z| < 0.9
clear n z; syms n z
assume(n,'integer');           % left-sided, no nonnegativity assumption

c = sym(9)/10;                 % exact 0.9
x_left = - (c^n) * heaviside(-n-1);
Xc = ztrans(x_left, n, z);
disp('A2(c) X(z) left-sided:'); pretty(Xc)
roc_c = '|z| < 0.9';

% Save ROC summary
fid = fopen(fullfile(outdir,'A2_ROC_summary.txt'),'w');
fprintf(fid,'A2(a): X(z)=1/(1-0.6 z^{-1}), ROC %s\n',roc_a);
fprintf(fid,'A2(b): X(z)=1/(1+0.8 z^{-1}), ROC %s\n',roc_b);
fprintf(fid,'A2(c): X(z)=-(z)/(z-0.9), ROC %s (left-sided)\n',roc_c);
fclose(fid);

%% ============== A3: Linearity & Time Shifting ==========================
% x1[n] = (0.5)^n u[n], x2[n] = (-0.5)^n u[n]
clear n z; syms n z
assume(n,'integer'); assumeAlso(n>=0);

x1 = (sym(1)/2)^n * heaviside(n);
x2 = (-sym(1)/2)^n * heaviside(n);

% (a) Z{ 2 x1[n] - 3 x2[n] }
X_lin = ztrans(2*x1 - 3*x2, n, z);
disp('A3(a) X_lin(z) = '); pretty(X_lin)

% (b) Time shift: Z{ x1[n-3] } = z^{-3} X1(z), ROC unchanged for right-sided
X_shift = ztrans(subs(x1, n, n-3), n, z);
disp('A3(b) X_shift(z) = '); pretty(X_shift)

% Save symbolic results
fid = fopen(fullfile(outdir,'A3_symbolic.txt'),'w');
fprintf(fid,'X_lin(z)   = %s\n', char(X_lin));
fprintf(fid,'Z{x1[n-3]} = %s\n', char(X_shift));
fclose(fid);

%% =================== A4: Inverse Z-Transform ===========================
clear n z; syms n z

% (a) X(z) = 1/(1 - 0.7 z^{-1})  -> x[n] = (0.7)^n u[n]
Xa = 1 / (1 - (sym(7)/10)*z^(-1));
xa = iztrans(Xa, z, n);
disp('A4(a) x[n] ='); pretty(xa)

% (b) X(z) = (1 - 0.5 z^{-1})/(1 - 0.8 z^{-1})
%     -> x[n] = (0.8)^n u[n] - 0.5 (0.8)^{n-1} u[n-1]
Xb = (1 - (sym(1)/2)*z^(-1)) / (1 - (sym(4)/5)*z^(-1));
xb = iztrans(Xb, z, n);
disp('A4(b) x[n] ='); pretty(xb)

fid = fopen(fullfile(outdir,'A4_iztrans.txt'),'w');
fprintf(fid,'xa[n] = %s\n', char(xa));
fprintf(fid,'xb[n] = %s\n', char(xb));
fclose(fid);

%% ===== A5: H(z), Poles/Zeros & Frequency Response (plots saved) ========
b = [1 -2.4  2.88];
a = [1 -0.8  0.64];

% (a) Pole-zero plot
f1 = figure('Name','A5_zplane');
zplane(b,a); grid on; title('A5: Pole–Zero Plot');
saveas(f1, fullfile(outdir,'A5_zplane.png'));

% Numeric poles/zeros/gain
zz = roots(b);    % zeros
pp = roots(a);    % poles
kk = b(1)/a(1);   % gain (here 1)

% Save values neatly
fid = fopen(fullfile(outdir,'A5_pz_values.txt'),'w');
fprintf(fid,'Zeros:\n');
for k = 1:numel(zz)
    fprintf(fid,'  % .6f %+.6fj (|r|=% .6f, angle=% .3f deg)\n', real(zz(k)), imag(zz(k)), abs(zz(k)), rad2deg(angle(zz(k))));
end
fprintf(fid,'Poles:\n');
for k = 1:numel(pp)
    fprintf(fid,'  % .6f %+.6fj (|r|=% .6f, angle=% .3f deg)\n', real(pp(k)), imag(pp(k)), abs(pp(k)), rad2deg(angle(pp(k))));
end
fprintf(fid,'Gain: %.6f\n', kk);
fclose(fid);

% (b) Magnitude & Phase response
[H,w] = freqz(b,a,512);  % 0..pi
f2 = figure('Name','A5_freqz');
subplot(2,1,1);
plot(w/pi, abs(H)); grid on; ylabel('|H(e^{j\omega})|'); title('Magnitude Response');
subplot(2,1,2);
plot(w/pi, angle(H)); grid on; xlabel('\omega/\pi'); ylabel('Phase (rad)');
saveas(f2, fullfile(outdir,'A5_freq_response.png'));

% (Optional) quick time-domain test signal
n = 0:511;
x = sin(0.2*pi*n) + 0.5*sin(0.8*pi*n);
y = filter(b,a,x);
f3 = figure('Name','A5_optional_time');
plot(n, x, 'DisplayName','x[n]'); hold on;
plot(n, y, 'DisplayName','y[n]');
grid on; legend; xlabel('n'); title('Optional Test: Input vs Output');
saveas(f3, fullfile(outdir,'A5_optional_time.png'));

disp('All done. Check the ./figures folder for outputs.');
