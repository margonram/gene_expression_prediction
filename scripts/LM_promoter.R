setwd("yourpath")

# Load packages
library(caret)
library(viridis)

# Upload and prepare data
data<-read.table("results/matrix_observations_promoters.txt",header = T) # substitute by your file
head(data,n=5)
data<-log2(data+0.1)
head(data,n=5)
data<-data.frame(data)
nrow(data)
ncol(data)

# Create sets of training and test
set.seed(47)
trainIndex <- createDataPartition(data$expression, p = .8,list = FALSE,times = 1)
head(trainIndex,n=2)
tail(trainIndex,n=2)
dataTrain <- data[ trainIndex,]
dataTest  <- data[-trainIndex,]
nrow(dataTrain)
nrow(dataTest)
head(dataTrain,n=2)
head(dataTest,n=2)

# Create random dataTrain
expression_rd <- sample(dataTrain$expression)
dataTrain_rd <- dataTrain
dataTrain_rd$expression <- expression_rd
head(dataTrain_rd,n=2)

# define training control: 10X cross-validation
train_control <- trainControl(method="repeatedcv", number=10, repeats = 3)

# learn models
model <- train(expression~H3K36me3+H3K4me3+H3K27me3+H3K27ac, # substitute by your formula
               data=dataTrain, trControl=train_control, method="lm")

model_rd <- train(expression~H3K36me3+H3K4me3+H3K27me3+H3K27ac, # substitute by your formula
                  data=dataTrain_rd, trControl=train_control, method="lm")

# save results
model
model_print<-capture.output(print(model))
out_print = paste("results/Promoter_print.txt", sep = "")
writeLines(model_print,out_print)
summary(model)
model_sum<-capture.output(summary(model))
out_sum = paste("results/Promoter_summary.txt", sep = "")
writeLines(model_sum,out_sum)

# calculate variable importance
varImp(model, scale = FALSE)
model_imp<-capture.output(varImp(model, scale = FALSE))
out_imp = paste("results/Promoter_varImp.txt", sep = "")
writeLines(model_imp,out_imp)

# calculate Pearson's correlation
exp_pred<-predict(model,dataTest)
r<-round(cor(exp_pred,dataTest$expression),2)
r
exp_pred_rd<-predict(model_rd,dataTest)
r_rd<-round(cor(exp_pred_rd,dataTest$expression),2)
r_rd

# plot predicted vs. measured
predicted<-c(exp_pred,exp_pred_rd)
measured<-c(dataTest$expression,dataTest$expression)
models<-c(rep("Model",nrow(dataTest)),rep("Random model",nrow(dataTest)))

text<-data.frame(label=c(paste("r = ",r),paste("r = ",r_rd)),models=c("Model","Random model"))
c<-data.frame(predicted,measured,models)

pdf("plots/Promoter_prediction.pdf", width = 10,height = 5)
ggplot(c) +
  geom_hex(aes(measured, predicted), bins = 100) +
  scale_fill_gradientn("", colours = rev(viridis(300)))+
  geom_smooth(aes(measured, predicted),method = "lm",level=0)+
  labs(title="Expression prediction in the test subset",x="Measured expression (log(FPKMs + 0.1))", y = "Predicted expression") +
  geom_text(data = text, mapping = aes(x = -Inf, y = -Inf, label = label),hjust = -2, vjust = -1, size=7) +
  theme_bw() +
  theme(legend.position="right",axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),legend.title=element_text(size=20),plot.title = element_text(size=24,face="bold"),
        strip.text = element_text(size=20,face="bold")) +
  facet_wrap(~ models)
dev.off()

# plot variable importance
out_varImp<-paste("plots/Promoter_varImp.pdf",sep="")
p<-ggplot(varImp(model,scale = FALSE))+
  geom_bar(stat="identity", fill="darkviolet")+
  labs(title="Promoter model variable importance",x="", y = "importance") +
  theme_bw() +
  theme(legend.position="bottom",axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"),
        plot.title = element_text(size=15,face="bold"),strip.text = element_text(size=12,face="bold")) 
ggsave(out_varImp,plot=p,device = "pdf",width = 5, height = 4)
