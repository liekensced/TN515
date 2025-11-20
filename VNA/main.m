clear all; clc

%% --- S11 (LightGrey) ---
M1 = readmatrix('S11LightGrey.csv', 'Delimiter',';', 'CommentStyle','!', 'DecimalSeparator',',');
f1 = M1(:,1)/1e6;
S11_dB = 20*log10(hypot(M1(:,2),M1(:,3)));
Gamma1 = 10.^(S11_dB/20);
VSWR1  = (1+Gamma1)./(1-Gamma1);

%% --- S22 (Brick) ---
M2 = readmatrix('S22Brick.csv', 'Delimiter',';', 'CommentStyle','!', 'DecimalSeparator',',');
f2 = M2(:,1)/1e6;
S22_dB = 20*log10(hypot(M2(:,2),M2(:,3)));
Gamma2 = 10.^(S22_dB/20);
VSWR2  = (1+Gamma2)./(1-Gamma2);

%% --- Plot both one above the other ---
figure('Position',[100 100 900 600]);

% Top: S11
subplot(2,1,1)
yyaxis left
plot(f1, S11_dB, '-o','LineWidth',1.5,'MarkerFaceColor','b')
ylabel('S11 (dB)'); ylim([-60 0])

yyaxis right
plot(f1, VSWR1, '-s','Color','#D95319','LineWidth',1.5,'MarkerFaceColor','#D95319')
ylabel('VSWR'); ylim([1 10])

grid on; title('S11 - LightGrey Antenna')
xlabel('Frequency (MHz)'); xlim([min(f1) max(f1)])

% Bottom: S22
subplot(2,1,2)
yyaxis left
plot(f2, S22_dB, '-o','LineWidth',1.5,'MarkerFaceColor','b')
ylabel('S22 (dB)'); ylim([-60 0])

yyaxis right
plot(f2, VSWR2, '-s','Color','#D95319','LineWidth',1.5,'MarkerFaceColor','#D95319')
ylabel('VSWR'); ylim([1 10])

grid on; title('S22 - Brick Antenna')
xlabel('Frequency (MHz)'); xlim([min(f2) max(f2)])

sgtitle('S11 vs S22 Comparison   (dual axis: dB left, VSWR right)')