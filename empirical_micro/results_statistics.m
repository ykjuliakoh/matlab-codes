%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RESULTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

true_beta=[0.5;1;1;0.5;1;1];
MC_N=500;
GII_LR_nonan=GII_LR_QN1(:,~all(isnan(GII_LR_QN1)));
[K,M_LR]=size(GII_LR_nonan);

MEAN_SML1=sum(SML_R30,2)/MC_N;
SD_SML1=std(SML_R30,0,2);
MSE_SML1=sum(((SML_R30-true_beta*ones(1,500)).^2),2)/MC_N;

MEAN_SML2=sum(SML_R50,2)/MC_N;
SD_SML2=std(SML_R50,0,2);
MSE_SML2=sum(((SML_R50-true_beta*ones(1,500)).^2),2)/MC_N;

MEAN_SML3=sum(SML_R100,2)/MC_N;
SD_SML3=std(SML_R100,0,2);
MSE_SML3=sum(((SML_R100-true_beta*ones(1,500)).^2),2)/MC_N;

MEAN_MSM1=sum(MSM_R30,2)/MC_N;
SD_MSM1=std(MSM_R30,0,2);
MSE_MSM1=sum(((MSM_R30-true_beta*ones(1,500)).^2),2)/MC_N;

MEAN_MSM2=sum(MSM_R50,2)/MC_N;
SD_MSM2=std(MSM_R50,0,2);
MSE_MSM2=sum(((MSM_R50-true_beta*ones(1,500)).^2),2)/MC_N;

MEAN_MSM3=sum(MSM_R100,2)/MC_N;
SD_MSM3=std(MSM_R100,0,2);
MSE_MSM3=sum(((MSM_R100-true_beta*ones(1,500)).^2),2)/MC_N;

MSE_GII_LR=sum(((GII_LR_nonan-true_beta*ones(1,M_LR)).^2),2)/M_LR;
MEAN_GII_LR=sum(GII_LR_nonan,2)/M_LR;
SD_GII_LR=std(GII_LR_nonan,0,2);

statistics_ALL=zeros(6,7,3);
statistics_ALL(:,:,1)=[MEAN_SML1,MEAN_SML2,MEAN_SML3,MEAN_MSM1,MEAN_MSM2,MEAN_MSM3,MEAN_GII_LR];
statistics_ALL(:,:,2)=[SD_SML1,SD_SML2,SD_SML3,SD_MSM1,SD_MSM2,SD_MSM3,SD_GII_LR];
statistics_ALL(:,:,3)=[MSE_SML1,MSE_SML2,MSE_SML3,MSE_MSM1,MSE_MSM2,MSE_MSM3,MSE_GII_LR];