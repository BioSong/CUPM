library(survival)
library(survminer)
load("./CUPM/results/V2_202507/4_CuScore_CMS/TCGA&ARGO_data.RData")
#cuScore1,cuScore2,comp_result11,comp_result21,tcga_cli,argo_cli
#-----------------------TCGA
df_tcga=cbind(tcga_cli[,c("sample","days_to_death","vital_status")],cuScore1$Final_Score)
colnames(df_tcga)[2:4]=c("time","status","score")
best_thres <- surv_cutpoint(df_tcga,
                                     time = "time",  # 生存时间列名
                                     event = "status",      # 生存事件列名
                                     variables = "score",  # 需要寻找阈值的变量列名
                                     minprop = 0.3,     # 最小比例，防止找到的阈值过于极端
                                     progressbar = TRUE)  # 显示进度条

df_tcga$group <- ifelse(df_tcga$score > best_thres$cutpoint$cutpoint, "High", "Low")
# 2. 拟合生存对象
surv_obj <- Surv(time = df_tcga$time, event = df_tcga$status)
# 3. 绘制KM曲线
fit <- survfit(surv_obj ~ group, data = df_tcga)
# 4. 可视化
ggsurvplot(fit,
           data = df_tcga,
           pval = TRUE, 
           risk.table = TRUE, 
           conf.int = TRUE, 
           palette = c("#E64B35", "#4DBBD5"),
           legend.title = "Cuproptosis Score",
           xlab = "Time (months)",
           ylab = "Survival Probability")
#---------------------------------------ARGO
df_argo=cbind(argo_cli[,c("sample","OS_days","Death")],cuScore2$Final_Score)
colnames(df_argo)[2:4]=c("time","status","score")
best_thres <- surv_cutpoint(df_argo,
                                     time = "time",  # 生存时间列名
                                     event = "status",      # 生存事件列名
                                     variables = "score",  # 需要寻找阈值的变量列名
                                     minprop = 0.3,     # 最小比例，防止找到的阈值过于极端
                                     progressbar = TRUE)  # 显示进度条


df_argo$group <- ifelse(df_argo$score > median(df_argo$score), "High", "Low")

# 2. 拟合生存对象
surv_obj <- Surv(time = df_argo$time, event = df_argo$status)

# 3. 绘制KM曲线
fit <- survfit(surv_obj ~ group, data = df_argo)

# 4. 可视化
ggsurvplot(fit,
           data = df_argo,
           pval = TRUE, 
           risk.table = TRUE, 
           conf.int = TRUE, 
           palette = c("#E64B35", "#4DBBD5"),
           legend.title = "Cuproptosis Score",
           xlab = "Time (months)",
           ylab = "Survival Probability")
#--------------------------combine tcga and argo
df_all=rbind(df_tcga,df_argo)
# 2. 拟合生存对象
surv_obj <- Surv(time = df_all$time, event = df_all$status)

# 3. 绘制KM曲线
fit <- survfit(surv_obj ~ group, data = df_all)

# 4. 可视化
ggsurvplot(fit,
           data = df_all,
           pval = TRUE, 
           risk.table = TRUE, 
           conf.int = TRUE, 
           palette = c("#E64B35", "#4DBBD5"),
           legend.title = "Cuproptosis Score",
           xlab = "Time (months)",
           ylab = "Survival Probability")

