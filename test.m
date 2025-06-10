%--------------Load dependency and image data----------------
load('RGB_absor.mat');
% singleH: stain absorption coefficient matrix
A = singleH';

load('Calibration.mat');
% f_E/f_H/f_L/f_O: calibration parameters

% Load RGB image
img = imread('Case05_E-L_04_RGB.jpg');
bgX2 = im2double(img);

load('Case05_E-L_04_GT.mat');
% GT: ground truth of stain abundance (scaled to 0-1)


%------------Calculate optical density from RGB image-------------
[H,W,band] = size(GT);
bgX1 = reshape(bgX2,W*H,3);
bgX1_s = round(sum(bgX1,2),2);
b99_t = tabulate(bgX1_s);
[b99_p,b99_in] = max(b99_t(:,2));
bgX1_s1 = b99_t(:,1);
b9_99ind = find( bgX1_s>bgX1_s1(b99_in) );

% Calculate RGB intensity of background
back_RGB = bgX1(b9_99ind,:);
back_RGB = mean(back_RGB);
back_R = back_RGB(1);
back_G = back_RGB(2);
back_B = back_RGB(3);

% Calculate optical density using Beer-Lambert Law
ODr(:,:,1) = -log10(bgX2(:,:,1)./back_R);
ODr(:,:,2) = -log10(bgX2(:,:,2)./back_G);
ODr(:,:,3) = -log10(bgX2(:,:,3)./back_B);
ODr(ODr<0) = 0;
ODx = reshape(ODr,W*H,3);
ODx(isinf(ODx)) = nan;
max_x = max(ODx);
ODx(isnan(ODx(:,1)),1) = max_x(1);
ODx(isnan(ODx(:,2)),2) = max_x(2);
ODx(isnan(ODx(:,3)),3) = max_x(3);
OD = reshape(ODx,H,W,3);

%--------------Select ROI---------------
ab = GT(431:880,501:1050,:);
od = OD(431:880,501:1050,:);
bgrgb = bgX2(431:880,501:1050,:);
[h,w,band] = size(ab);
odx = reshape(od,h*w,3);
ab1 = reshape(ab,h*w,4);
ab2 = ab1';

% RGB image
figure
imshow(bgrgb) 

stain_names = {'EY','H','LG','OG'};
% Stain abundance map of ground truth
figure
set(gcf,'Position',[0,0,w+100,h]);
for i = 1:4
    subplot(2,2,i)
    hm1 = heatmap(ab(:,:,i),'GridVisible','off','CellLabelColor','none');
    hm1.Colormap = jet;
    clim([0 1])
    cd1 = hm1.XDisplayLabels;                                    % Current Display Labels
    hm1.XDisplayLabels = repmat(' ',size(cd1,1), size(cd1,2));   % Blank Display Labels
    cd2 = hm1.YDisplayLabels;                                    % Current Display Labels
    hm1.YDisplayLabels = repmat(' ',size(cd2,1), size(cd2,2));   % Blank Display Labels
    title({stain_names{i};' '})
end

%% Stain Unmixing with nonnegativity, Weight Nucleus Sparsity and Total Variation

Y = odx';

% Parameters (Can be modified)
gamma = 0.01;
maxIter = 2000;
Lam1 = 2e-6;
LamTV = 1e-3;

[X,Res_iX] = SUnWNS_TV(A,Y,'MU',gamma,'POSITIVITY','yes', ...
                  'LAMBDA_1',Lam1,'LAMBDA_TV', LamTV, ...
                  'IM_SIZE',[h,w],'AL_ITERS',maxIter, 'VERBOSE','yes');

Xp = X';
Xp = reshape(Xp, [h,w,4]);
Xp(Xp<0) = 0;
X1 = reshape(Xp,h*w,4);
X2 = X1;

% Calibration
X2(:,1) = f_E(X2(:,1));
X2(:,2) = f_H(X2(:,2));
X2(:,3) = f_L(X2(:,3));
X2(:,4) = f_O(X2(:,4));
XD = reshape(X2,h,w,4);

% Calculate SRE and RMSE
P_XX = sum(sum((ab1-X2).^2));
SRE_XX = 10*log10(Pab/P_XX);
P_X = sum(sum((ab1-X2).^2));
SRE_X = 10*log10(Pab/P_X);
RMSE_X = rmse(X2,ab1,"all");
PS_E_E = sum(sum((ab1(:,1)-X2(:,1)).^2));
PS_E_H = sum(sum((ab1(:,2)-X2(:,2)).^2));
PS_E_L = sum(sum((ab1(:,3)-X2(:,3)).^2));
PS_E_O = sum(sum((ab1(:,4)-X2(:,4)).^2));
SRE_X4(1) = 10*log10(P_E/PS_E_E);
SRE_X4(2) = 10*log10(P_H/PS_E_H);
SRE_X4(3) = 10*log10(P_L/PS_E_L);
SRE_X4(4) = 10*log10(P_O/PS_E_O);
RMSE_X4(1) = rmse(X2(:,1),ab1(:,1));
RMSE_X4(2) = rmse(X2(:,2),ab1(:,2));
RMSE_X4(3) = rmse(X2(:,3),ab1(:,3));
RMSE_X4(4) = rmse(X2(:,4),ab1(:,4));

% Estimated stain abundance map
figure
set(gcf,'Position',[0,0,w,h]);
for i = 1:4
    subplot(2,2,i)
    hm1 = heatmap(XD(:,:,i),'GridVisible','off','CellLabelColor','none');
    hm1.Colormap = jet;
    clim([0 1])
    cd1 = hm1.XDisplayLabels;                                    % Current Display Labels
    hm1.XDisplayLabels = repmat(' ',size(cd1,1), size(cd1,2));   % Blank Display Labels
    cd2 = hm1.YDisplayLabels;                                    % Current Display Labels
    hm1.YDisplayLabels = repmat(' ',size(cd2,1), size(cd2,2));   % Blank Display Labels
    title({stain_names{i};['SRE=',num2str(SRE_X4(i),'%.4f'),'  RMSE=',num2str(RMSE_X4(i),'%.4f')]})
end
sgtitle({'SUnWNS';['SRE=',num2str(SRE_X,'%.4f'),'  RMSE=',num2str(RMSE_X,'%.4f')]})



% %% Traditional unmixing 
% %  Color Deconvolution using Moore-Penrose pseudo-inverse
% H0 = singleH;
% Hinv1 = pinv(H0);
% C2_1 = odx*Hinv1;
% 
% % Calibration
% C2_1(:,1) = f_E(C2_1(:,1));
% C2_1(:,2) = f_H(C2_1(:,2));
% C2_1(:,3) = f_L(C2_1(:,3));
% C2_1(:,4) = f_O(C2_1(:,4));
% C2_1(C2_1<0)=0.0;
% C1 = reshape(C2_1,[h,w,4]);
% 
% % Calculate SRE and RMSE
% Pab = sum(sum(ab1.^2));
% P_C = sum(sum((ab1-C2_1).^2));
% SRE_C = 10*log10(Pab/P_C);
% RMSE_C = rmse(C2_1,ab1,"all");
% P_E = sum(sum((ab1(:,1)).^2));
% P_H = sum(sum((ab1(:,2)).^2));
% P_L = sum(sum((ab1(:,3)).^2));
% P_O = sum(sum((ab1(:,4)).^2));
% PC_E_E = sum(sum((ab1(:,1)-C2_1(:,1)).^2));
% PC_E_H = sum(sum((ab1(:,2)-C2_1(:,2)).^2));
% PC_E_L = sum(sum((ab1(:,3)-C2_1(:,3)).^2));
% PC_E_O = sum(sum((ab1(:,4)-C2_1(:,4)).^2));
% SRE_C4(1) = 10*log10(P_E/PC_E_E);
% SRE_C4(2) = 10*log10(P_H/PC_E_H);
% SRE_C4(3) = 10*log10(P_L/PC_E_L);
% SRE_C4(4) = 10*log10(P_O/PC_E_O);
% RMSE_C4(1) = rmse(C2_1(:,1),ab1(:,1));
% RMSE_C4(2) = rmse(C2_1(:,2),ab1(:,2));
% RMSE_C4(3) = rmse(C2_1(:,3),ab1(:,3));
% RMSE_C4(4) = rmse(C2_1(:,4),ab1(:,4));
% 
% % Estimated stain abundance map
% figure
% set(gcf,'Position',[0,0,w,h]);
% for i = 1:4
%     subplot(2,2,i)
%     hm1 = heatmap(C1(:,:,i),'GridVisible','off','CellLabelColor','none');
%     hm1.Colormap = jet;
%     clim([0 1])
%     cd1 = hm1.XDisplayLabels;                                    % Current Display Labels
%     hm1.XDisplayLabels = repmat(' ',size(cd1,1), size(cd1,2));   % Blank Display Labels
%     cd2 = hm1.YDisplayLabels;                                    % Current Display Labels
%     hm1.YDisplayLabels = repmat(' ',size(cd2,1), size(cd2,2));   % Blank Display Labels
%     title({stain_names{i};['SRE=',num2str(SRE_C4(i),'%.4f'),'  RMSE=',num2str(RMSE_C4(i),'%.4f')]})
% end
% sgtitle({'Color Deconvolution';['SRE=',num2str(SRE_C,'%.4f'),'  RMSE=',num2str(RMSE_C,'%.4f')]})
