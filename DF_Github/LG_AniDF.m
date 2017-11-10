%%
function [AniDF]=LG_AniDF(DF_out,sample,isCreateAVI)
N = length(DF_out.xP);
AniDF.xP = zeros(N,1);
AniDF.yP = zeros(N,1);
AniDF.zP = zeros(N,1);
AniDF.dxP = zeros(N,1);
AniDF.dyP = zeros(N,1);
AniDF.dzP = zeros(N,1);

AniDF.xQA = zeros(N,1);
AniDF.yQA = zeros(N,1);
AniDF.zQA = zeros(N,1);
AniDF.xQB = zeros(N,1);
AniDF.yQB = zeros(N,1);
AniDF.zQB = zeros(N,1);

AniDF.phiA = zeros(N,1);
AniDF.thetaA = zeros(N,1);
AniDF.psiA = zeros(N,1);
AniDF.phiB = zeros(N,1);
AniDF.thetaB = zeros(N,1);
AniDF.psiB = zeros(N,1);

AniDF.dphiA = zeros(N,1);
AniDF.dthetaA = zeros(N,1);
AniDF.dpsiA = zeros(N,1);
AniDF.dphiB = zeros(N,1);
AniDF.dthetaB = zeros(N,1);
AniDF.dpsiB = zeros(N,1);

AniDF.alpha1 = zeros(N,1);
AniDF.beta1 = zeros(N,1);
AniDF.alpha2 = zeros(N,1);
AniDF.beta2 = zeros(N,1);

AniDF.dalpha1 = zeros(N,1);
AniDF.dbeta1 = zeros(N,1);
AniDF.dalpha2 = zeros(N,1);
AniDF.dbeta2 = zeros(N,1);


AniDF.fA = zeros(N,1);
AniDF.Mx1 = zeros(N,1);
AniDF.My1 = zeros(N,1);
AniDF.Mz1 = zeros(N,1);
AniDF.fB = zeros(N,1);
AniDF.Mx2 = zeros(N,1);
AniDF.My2 = zeros(N,1);
AniDF.Mz2 = zeros(N,1);

Jac = @(phi,theta)([1,0,-sin(theta);...
    0,cos(phi),sin(phi)*cos(theta);...
    0,-sin(phi),cos(phi)*cos(theta)]);

for i = 1:N
    %% alpha1,beat1
    AniDF.alpha1(i) =  (DF_out.aA0(1,i)); % Î´¸üÐÂ
    AniDF.dalpha1(i) = (DF_out.aA1(1,i));
    
    AniDF.beta1(i) = (DF_out.bA0(1,i));
    AniDF.dbeta1(i) = (DF_out.bA1(1,i));
    %% alpha2,beta2
    AniDF.alpha2(i) =  (DF_out.aB0(1,i));
    AniDF.beta2(i) =   (DF_out.bB0(1,i));
    AniDF.dalpha2(i) = (DF_out.aB1(1,i));
    AniDF.dbeta2(i) =  (DF_out.bB1(1,i));
    %%
    AniDF.xQA(i) = DF_out.xQ1(i);
    AniDF.yQA(i) = DF_out.yQ1(i);
    AniDF.zQA(i) = DF_out.zQ1(i);
    
    AniDF.xQB(i) = DF_out.xQ2(i);
    AniDF.yQB(i) = DF_out.yQ2(i);
    AniDF.zQB(i) = DF_out.zQ2(i);
    %%
    AniDF.xP(i) = DF_out.xP(i);
    AniDF.yP(i) = DF_out.yP(i);
    AniDF.zP(i) = DF_out.zP(i);
    
    AniDF.dxP(i) = DF_out.dxP(i);
    AniDF.dyP(i) = DF_out.dyP(i);
    AniDF.dzP(i) = DF_out.dzP(i);
    %%
    AniDF.phiA(i) = DF_out.eulerA(1,i);
    AniDF.thetaA(i) = DF_out.eulerA(2,i);
    AniDF.psiA(i) = DF_out.eulerA(3,i);
    
    temp = (Jac(AniDF.phiA(i),AniDF.thetaA(i)))\DF_out.w1(:,i);
    AniDF.dphiA(i) = temp(1);
    AniDF.dthetaA(i) = temp(2);
    AniDF.dpsiA(i) = temp(3);
    
    
    AniDF.phiB(i) = DF_out.eulerB(1,i);
    AniDF.thetaB(i) = DF_out.eulerB(2,i);
    AniDF.psiB(i) = DF_out.eulerB(3,i);
    
    temp = Jac(AniDF.phiB(i),AniDF.thetaB(i))\DF_out.w2(:,i);
    AniDF.dphiB(i) = temp(1);
    AniDF.dthetaB(i) = temp(2);
    AniDF.dpsiB(i) = temp(3);    
    %%
    AniDF.fA(i) = DF_out.fA(i);
    AniDF.fB(i) = DF_out.fB(i);    
    %%
    AniDF.Mx1(i) = DF_out.MA(1,i);
    AniDF.My1(i) = DF_out.MA(2,i);
    AniDF.Mz1(i) = DF_out.MA(3,i);
    %
    AniDF.Mx2(i) = DF_out.MB(1,i);
    AniDF.My2(i) = DF_out.MB(2,i);
    AniDF.Mz2(i) = DF_out.MB(3,i);
end


posQ1 = [AniDF.xQA,AniDF.yQA,AniDF.zQA];
posQ2 = [AniDF.xQB,AniDF.yQB,AniDF.zQB];
posP =  [AniDF.xP,AniDF.yP,AniDF.zP];

R1 = euler2rotMat(AniDF.phiA,AniDF.thetaA,AniDF.psiA);
R2 = euler2rotMat(AniDF.phiB,AniDF.thetaB,AniDF.psiB);

SamplePlotFreq = sample;
% SamplePlotFreq = 25*2000;
% isCreateAVI = 1;
isFixView = false;
LG_6DoFAnimation_door(posQ1,R1,posQ2,R2,posP,SamplePlotFreq,'DotsOnly',isCreateAVI,isFixView);
end
% hold on;
% plot3(des_x,des_y,0*des_z,'r');