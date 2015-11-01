
library("ggplot2")
library("reshape2")

summary(ts_record_2randelecs_all)
summary(ts_record_4randelecs_all)
summary(ts_record_6randelecs_all)
summary(ts_record_8randelecs_all)

ts_record_2randelecs_all$chans = "two"
ts_record_4randelecs_all$chans = "four"
ts_record_6randelecs_all$chans = "six"
ts_record_8randelecs_all$chans = "eight"

ts_record_allchanperms = rbind(ts_record_2randelecs_all,ts_record_4randelecs_all,ts_record_6randelecs_all,ts_record_8randelecs_all)
ts_record_allchanperms$chans <- factor(ts_record_allchanperms$chans, levels = c("two","four","six","eight"))
summary(ts_record_allchanperms)

ggplot(ts_record_allchanperms, aes(x=randelec, fill=chans)) + geom_density(alpha=0.3)+coord_cartesian(xlim=c(-5,5))
ggplot(ts_record_allchanperms, aes(x=fixedelec, fill=chans)) + geom_density(alpha=0.3)+coord_cartesian(xlim=c(-5,5))
ggplot(ts_record_allchanperms, aes(x=avgelec, fill=chans)) + geom_density(alpha=0.3)+coord_cartesian(xlim=c(-5,5))

sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="two",]$randelec>2.048 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="two",]$randelec < (-2.048))/100 #false positive rate for random electrode, using Satterthwaite's approximation to degrees of freedom
sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="four",]$randelec>2.048 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="four",]$randelec < (-2.048))/100 #false positive rate for random electrode, using Satterthwaite's approximation to degrees of freedom
sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="six",]$randelec>2.048 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="six",]$randelec < (-2.048))/100 #false positive rate for random electrode, using Satterthwaite's approximation to degrees of freedom
sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="eight",]$randelec>2.048 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="eight",]$randelec < (-2.048))/100 #false positive rate for random electrode, using Satterthwaite's approximation to degrees of freedom

sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="two",]$fixedelec>1.993 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="two",]$fixedelec < (-1.993))/100 #false positive rate for fixed electrode
sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="four",]$fixedelec>1.993 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="four",]$fixedelec < (-1.993))/100 #false positive rate for fixed electrode
sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="six",]$fixedelec>1.993 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="six",]$fixedelec < (-1.993))/100 #false positive rate for fixed electrode
sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="eight",]$fixedelec>1.993 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="eight",]$fixedelec < (-1.993))/100 #false positive rate for fixed electrode

sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="two",]$avgelec>2.060 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="two",]$avgelec < (-2.060))/100 #false positive rate for averaged electrode
sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="four",]$avgelec>2.060 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="four",]$avgelec < (-2.060))/100 #false positive rate for averaged electrode
sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="six",]$avgelec>2.060 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="six",]$avgelec < (-2.060))/100 #false positive rate for averaged electrode
sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="eight",]$avgelec>2.060 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="eight",]$avgelec < (-2.060))/100 #false positive rate for averaged electrode


fprates_df <- data.frame(fprate= numeric(0), chans= character(0), model= character(0))
fprates_df <- rbind(fprates_df,list(fprate=I(sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="two",]$randelec>2.048 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="two",]$randelec < (-2.048))/100),chans="two",model="randelec"))
fprates_df$chans <- factor(fprates_df$chans, levels = c("two","four","six","eight"))
fprates_df$model <- factor(fprates_df$model, levels = c("randelec","fixedelec","avgelec"))
fprates_df <- rbind(fprates_df,list(fprate=I(sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="four",]$randelec>2.048 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="four",]$randelec < (-2.048))/100),chans="four",model="randelec"))
fprates_df <- rbind(fprates_df,list(fprate=I(sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="six",]$randelec>2.048 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="six",]$randelec < (-2.048))/100),chans="six",model="randelec"))
fprates_df <- rbind(fprates_df,list(fprate=I(sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="eight",]$randelec>2.048 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="eight",]$randelec < (-2.048))/100),chans="eight",model="randelec"))
fprates_df <- rbind(fprates_df,list(fprate=I(sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="two",]$fixedelec>1.993 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="two",]$fixedelec < (-1.993))/100),chans="two",model="fixedelec"))
fprates_df <- rbind(fprates_df,list(fprate=I(sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="four",]$fixedelec>1.993 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="four",]$fixedelec < (-1.993))/100),chans="four",model="fixedelec"))
fprates_df <- rbind(fprates_df,list(fprate=I(sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="six",]$fixedelec>1.993 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="six",]$fixedelec < (-1.993))/100),chans="six",model="fixedelec"))
fprates_df <- rbind(fprates_df,list(fprate=I(sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="eight",]$fixedelec>1.993 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="eight",]$fixedelec < (-1.993))/100),chans="eight",model="fixedelec"))
fprates_df <- rbind(fprates_df,list(fprate=I(sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="two",]$avgelec>2.060 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="two",]$avgelec < (-2.060))/100),chans="two",model="avgelec"))
fprates_df <- rbind(fprates_df,list(fprate=I(sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="four",]$avgelec>2.060 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="four",]$avgelec < (-2.060))/100),chans="four",model="avgelec"))
fprates_df <- rbind(fprates_df,list(fprate=I(sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="six",]$avgelec>2.060 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="six",]$avgelec < (-2.060))/100),chans="six",model="avgelec"))
fprates_df <- rbind(fprates_df,list(fprate=I(sum(ts_record_allchanperms[ts_record_allchanperms["chans"]=="eight",]$avgelec>2.060 | ts_record_allchanperms[ts_record_allchanperms["chans"]=="eight",]$avgelec < (-2.060))/100),chans="eight",model="avgelec"))

ggplot(fprates_df, aes(x = chans, y = fprate, color=model, group = model, label=fprate)) + 
  geom_point(size = 5)  + geom_text(hjust=0.75, vjust=-1) + geom_line(size = 1.25) + geom_line(y=5,color="black") + annotate("text",x="six",y=5.25,label="5% false positive")
