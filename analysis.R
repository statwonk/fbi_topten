

## Plotting

# define custom function to create a survival data.frame
createSurvivalFrame <- function(f.survfit){
  # initialise frame variable
  f.frame <- NULL
  # check if more then one strata
  if(length(names(f.survfit$strata)) == 0){
    # create data.frame with data from survfit
    f.frame <- data.frame(time=f.survfit$time, n.risk=f.survfit$n.risk, n.event=f.survfit$n.event, n.censor = f.survfit
                          $n.censor, surv=f.survfit$surv, upper=f.survfit$upper, lower=f.survfit$lower)
    # create first two rows (start at 1)
    f.start <- data.frame(time=c(0, f.frame$time[1]), n.risk=c(f.survfit$n, f.survfit$n), n.event=c(0,0),
                          n.censor=c(0,0), surv=c(1,1), upper=c(1,1), lower=c(1,1))
    # add first row to dataset
    f.frame <- rbind(f.start, f.frame)
    # remove temporary data
    rm(f.start)
  }
  else {
    # create vector for strata identification
    f.strata <- NULL
    for(f.i in 1:length(f.survfit$strata)){
      # add vector for one strata according to number of rows of strata
      f.strata <- c(f.strata, rep(names(f.survfit$strata)[f.i], f.survfit$strata[f.i]))
    }
    # create data.frame with data from survfit (create column for strata)
    f.frame <- data.frame(time=f.survfit$time, n.risk=f.survfit$n.risk, n.event=f.survfit$n.event, n.censor = f.survfit
                          $n.censor, surv=f.survfit$surv, upper=f.survfit$upper, lower=f.survfit$lower, strata=factor(f.strata))
    # remove temporary data
    rm(f.strata)
    # create first two rows (start at 1) for each strata
    for(f.i in 1:length(f.survfit$strata)){
      # take only subset for this strata from data
      f.subset <- subset(f.frame, strata==names(f.survfit$strata)[f.i])
      # create first two rows (time: 0, time of first event)
      f.start <- data.frame(time=c(0, f.subset$time[1]), n.risk=rep(f.survfit[f.i]$n, 2), n.event=c(0,0),
                            n.censor=c(0,0), surv=c(1,1), upper=c(1,1), lower=c(1,1), strata=rep(names(f.survfit$strata)[f.i],
                                                                                                 2))
      # add first two rows to dataset
      f.frame <- rbind(f.start, f.frame)
      # remove temporary data
      rm(f.start, f.subset)
    }
    # reorder data
    f.frame <- f.frame[order(f.frame$strata, f.frame$time), ]
    # rename row.names
    rownames(f.frame) <- NULL
  }
  # return frame
  return(f.frame)
}

# define custom function to draw kaplan-meier curve with ggplot
qplot_survival <- function(f.frame, f.CI="default", f.shape=3){
  # use different plotting commands dependig whether or not strata's are given
  if("strata" %in% names(f.frame) == FALSE){
    # confidence intervals are drawn if not specified otherwise
    if(f.CI=="default" | f.CI==TRUE ){
      # create plot with 4 layers (first 3 layers only events, last layer only censored)
      # hint: censoring data for multiple censoring events at timepoint are overplotted
      # (unlike in plot.survfit in survival package)
      ggplot(data=f.frame) + geom_step(aes(x=time, y=surv), direction="hv") + geom_step(aes(x=time,
                                                                                            y=upper), directions="hv", linetype=2) + geom_step(aes(x=time,y=lower), direction="hv", linetype=2) +
        geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)
    }
    else {
      # create plot without confidence intervalls
      ggplot(data=f.frame) + geom_step(aes(x=time, y=surv), direction="hv") +
        geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)
    }
  }
  else {
    if(f.CI=="default" | f.CI==FALSE){
      # without CI
      ggplot(data=f.frame, aes(group=strata, colour=strata)) + geom_step(aes(x=time, y=surv),
                                                                         direction="hv") + geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)
    }
    else {
      # with CI (hint: use alpha for CI)
      ggplot(data=f.frame, aes(colour=strata, group=strata)) + geom_step(aes(x=time, y=surv),
                                                                         direction="hv") + geom_step(aes(x=time, y=upper), directions="hv", linetype=2, alpha=0.5) +
        geom_step(aes(x=time,y=lower), direction="hv", linetype=2, alpha=0.5) +
        geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)
    }
  }
}

ggkm <- function(sfit,
                 table = TRUE,
                 returns = FALSE,
                 xlabs = "Time",
                 ylabs = "Survival Probability",
                 xlims = c(0,max(sfit$time)),
                 ylims = c(0,1),
                 ystratalabs = NULL,
                 ystrataname = NULL,
                 timeby = 100,
                 main = "Kaplan-Meier Plot",
                 pval = TRUE,
                 subs = NULL,
                 ...) {
  
  #############
  # libraries #
  #############
  
  #Check if the following packages have been installed. If not, install them
  if (!"ggplot2" %in% installed.packages()) install.packages("ggplot2")
  if (!"survival" %in% installed.packages()) install.packages("survival")
  if (!"gridExtra" %in% installed.packages()) install.packages("gridExtra")
  if (!"reshape" %in% installed.packages()) install.packages("reshape")
  
  suppressPackageStartupMessages(library(ggplot2, warn.conflicts=FALSE))
  suppressPackageStartupMessages(library(survival, warn.conflicts=FALSE))
  suppressPackageStartupMessages(library(gridExtra, warn.conflicts=FALSE))
  suppressPackageStartupMessages(library(reshape, warn.conflicts=FALSE))
  
  #################################
  # sorting the use of subsetting #
  #################################
  
  times <- seq(0, max(sfit$time), by = timeby)
  
  if(is.null(subs)){
    if(length(levels(summary(sfit)$strata)) == 0) {
      subs1 <- 1
      subs2 <- 1:length(summary(sfit,censored=T)$time)
      subs3 <- 1:length(summary(sfit,times = times,extend = TRUE)$time)
    } else {
      subs1 <- 1:length(levels(summary(sfit)$strata))
      subs2 <- 1:length(summary(sfit,censored=T)$strata)
      subs3 <- 1:length(summary(sfit,times = times,extend = TRUE)$strata)
    }
  } else{
    for(i in 1:length(subs)){
      if(i==1){
        ssvar <- paste("(?=.*\\b=",subs[i],sep="")
      }
      if(i==length(subs)){
        ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],"\\b)",sep="")
      }
      if(!i %in% c(1, length(subs))){
        ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],sep="")
      }
      if(i==1 & i==length(subs)){
        ssvar <- paste("(?=.*\\b=",subs[i],"\\b)",sep="")
      }
    }
    subs1 <- which(regexpr(ssvar,levels(summary(sfit)$strata), perl=T)!=-1)
    subs2 <- which(regexpr(ssvar,summary(sfit,censored=T)$strata, perl=T)!=-1)
    subs3 <- which(regexpr(ssvar,summary(sfit,times = times,extend = TRUE)$strata, perl=T)!=-1)
  }
  
  if( !is.null(subs) ) pval <- FALSE
  
  ##################################
  # data manipulation pre-plotting #
  ##################################
  
  if(length(levels(summary(sfit)$strata)) == 0) {
    #[subs1]
    if(is.null(ystratalabs)) ystratalabs <- as.character(sub("group=*","","All"))
  } else {
    #[subs1]
    if(is.null(ystratalabs)) ystratalabs <- as.character(sub("group=*","",names(sfit$strata)))
  }
  
  if(is.null(ystrataname)) ystrataname <- "Strata"
  m <- max(nchar(ystratalabs))
  times <- seq(0, max(sfit$time), by = timeby)
  
  if(length(levels(summary(sfit)$strata)) == 0) {
    Factor <- factor(rep("All",length(subs2)))
  } else {
    Factor <- factor(summary(sfit, censored = T)$strata[subs2])
  }
  
  #Data to be used in the survival plot
  .df <- data.frame(
    time = sfit$time[subs2],
    n.risk = sfit$n.risk[subs2],
    n.event = sfit$n.event[subs2],
    surv = sfit$surv[subs2],
    strata = Factor,
    upper = sfit$upper[subs2],
    lower = sfit$lower[subs2]
  )
  
  #Final changes to data for survival plot
  levels(.df$strata) <- ystratalabs
  zeros <- data.frame(time = 0, surv = 1,
                      strata = factor(ystratalabs, levels=levels(.df$strata)),
                      upper = 1, lower = 1)
  .df <- rbind.fill(zeros, .df)
  d <- length(levels(.df$strata))
  
  ###################################
  # specifying plot parameteres etc #
  ###################################
  
  p <- ggplot( .df, aes(time, surv)) +
    geom_step(aes(linetype = strata), size = 0.7) +
    theme_bw() +
    theme(axis.title.x = element_text(vjust = 0.5)) +
    scale_x_continuous(xlabs, breaks = times, limits = xlims) +
    scale_y_continuous(ylabs, limits = ylims) +
    theme(panel.grid.minor = element_blank()) +
    # MOVE LEGEND HERE BELOW [first is x dim, second is y dim]
    theme(legend.position = c(ifelse(m < 10, .85, .75),ifelse(d < 4, .85, .8))) +
    theme(legend.key = element_rect(colour = NA)) +
    labs(linetype = ystrataname) +
    theme(plot.margin = unit(c(0, 1, .5,ifelse(m < 10, 1.5, 2.5)),"lines")) +
    ggtitle(main)
  
  ## Create a blank plot for place-holding
  blank.pic <- ggplot(.df, aes(time, surv)) +
    geom_blank() + theme_bw() +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),panel.border = element_blank())
  
  #####################
  # p-value placement #
  #####################a
  
  if(length(levels(summary(sfit)$strata)) == 0) pval <- FALSE
  
  if(pval) {
    sdiff <- survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
    pval <- pchisq(sdiff$chisq,length(sdiff$n) - 1,lower.tail = FALSE)
    pvaltxt <- ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))
    # MOVE P-VALUE LEGEND HERE BELOW [set x and y]
    p <- p + annotate("text",x = 150, y = 0.1,label = pvaltxt)
  }
  
  ###################################################
  # Create table graphic to include at-risk numbers #
  ###################################################
  
  if(length(levels(summary(sfit)$strata)) == 0) {
    Factor <- factor(rep("All",length(subs3)))
  } else {
    Factor <- factor(summary(sfit,times = times,extend = TRUE)$strata[subs3])
  }
  
  if(table) {
    risk.data <- data.frame(
      strata = Factor,
      time = summary(sfit,times = times,extend = TRUE)$time[subs3],
      n.risk = summary(sfit,times = times,extend = TRUE)$n.risk[subs3]
    )
    risk.data$strata <- factor(risk.data$strata, levels=rev(levels(risk.data$strata)))
    
    data.table <- ggplot(risk.data,aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
      geom_text(size = 3.5) + theme_bw() +
      scale_y_discrete(breaks = as.character(levels(risk.data$strata)),
                       labels = rev(ystratalabs)) +
      scale_x_continuous("Numbers at risk", limits = xlims) +
      theme(axis.title.x = element_text(size = 10, vjust = 1),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),axis.text.x = element_blank(),
            axis.ticks = element_blank(),axis.text.y = element_text(face = "bold",hjust = 1))
    
    data.table <- data.table +
      theme(legend.position = "none") + xlab(NULL) + ylab(NULL)
    
    # ADJUST POSITION OF TABLE FOR AT RISK
    data.table <- data.table +
      theme(plot.margin = unit(c(-1.5, 1, 0.1, ifelse(m < 10, 2.5, 3.5) - 0.15 * m), "lines"))
    
    #######################
    # Plotting the graphs #
    #######################
    
    grid.arrange(p, blank.pic, data.table, clip = FALSE, nrow = 3,
                 ncol = 1, heights = unit(c(2, .1, .25),c("null", "null", "null")))
    
    if(returns) {
      a <- arrangeGrob(p, blank.pic, data.table, clip = FALSE, nrow = 3,
                       ncol = 1, heights = unit(c(2, .1, .25), c("null", "null", "null")))
      return(a)
    }#if
  } else {
    if(returns) return(p)
  }#else
}#function ggkm

df <- read.csv("data.csv",
                 header = T,
                 stringsAsFactors = F)

df$date_put_on_list <- as.POSIXct(df$date_put_on_list, format = "%m/%d/%Y", tz = "UTC")
df$follow_up_date <- as.POSIXct(df$follow_up_date, format = "%m/%d/%Y", tz = "UTC")
df$days_to_capture <- ifelse(is.na(df$follow_up_date) | df$follow_up_date == "",
                            difftime(Sys.time(), df$date_put_on_list, units = "days"),
                            difftime(df$follow_up_date, df$date_put_on_list, units = "days"))
df$censor <- ifelse(!is.na(df$follow_up_date), 1, 0)
mean(df$days_to_capture[df$days_to_capture > 0 & df$censor == 1], na.rm = T)


## Start #####
[1] 402.2239 # days

bootstrap_samples <- 10000
sample_size <- length(df$days_to_capture[df$days_to_capture > 0 & df$censor == 1]) * 0.8
means <- rep(NA, bootstrap_samples)
for(i in 1:bootstrap_samples) {
  only_fugitives_not_caught_before_annoucement <- df[df$days_to_capture > 0 & df$censor == 1, ]
  sample <- sample(1:length(only_fugitives_not_caught_before_annoucement$days_to_capture),
                   sample_size, replace = T)
  means[i] <- mean(only_fugitives_not_caught_before_annoucement$days_to_capture[sample], na.rm = T)
}

# install.packages("ggplot2"); library(ggplot2)
# install.packages("ggthemes"); library(ggthemes)
ggplot(as.data.frame(means), aes(x = means)) +
  geom_histogram(binwidth = 10, colour = "black", fill = "grey") +
  geom_vline(xintercept = quantile(as.data.frame(means)$means, 0.05), colour = "red") + # 5% quantile
  geom_vline(xintercept = quantile(as.data.frame(means)$means, 0.95), colour = "red") + # 95% quantile
  scale_x_continuous(breaks = seq(0, 10000, 100)) +
  theme_few() +
  ggtitle("Histogram of FBI Ten Most Wanted List Freedom") +
  ylab("Frequency of means by 10 day bin") +
  xlab("Average time free in days") +
  theme(axis.text = element_text(size = 16),
        plot.title = element_text(size = 18))

# difftime(Sys.time(), min(df$date_put_on_list[df$censor == 0]), units = "days")
# Time difference of 15892.45 days
# df$days_to_capture <- ifelse(df$id == 313, 15892.45, df$days_to_capture)
# mean(df$days_to_capture[df$days_to_capture > 0 & (df$censor == 1 | df$id == 313)], na.rm = T)

# df$year_put_on_list <- floor(as.integer(format(df$date_put_on_list, "%y"))/ 10)

df <- df[!is.na(df$censor) & df$days_to_capture >= 0, ]
# install.packages("survival"); library(survival)
print(fit <- survfit(Surv(days_to_capture, censor) ~ 1, data = df))
print(fit, rmean = "common")

plot(fit, fun = "cumhaz",
    # xmax = 365.25*3,
     ylab = "Rate of one fugitive being caught",
     xlab = "Days from being listed")
fit <- createSurvivalFrame(fit)

library(ggthemes)
few_mod <- function (base_size = 12, base_family = "") 
{
  colors <- ggthemes_data$few
  gray <- colors$medium["gray"]
  black <- colors$dark["black"]
  theme_bw(base_size = base_size, base_family = base_family) + 
    theme(line = element_line(colour = gray), rect = element_rect(fill = "white", 
                                                                  colour = NA), text = element_text(colour = black), 
          axis.ticks = element_line(colour = gray), legend.key = element_rect(colour = NA), 
          panel.border = element_rect(colour = gray), 
          strip.background = element_rect(fill = "white", colour = NA))
}

# install.packages("scales"); library(scales)
ggplot(fit[fit$time >= 0, ]) + 
  geom_step(aes(x = time / 365.25, y = surv), direction="hv", size = 1) +
  geom_ribbon(aes(x = time / 365.25, ymin = lower, ymax = upper), alpha = 0.2, colour = "blue") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 5)) +
  # scale_x_log10(breaks = c(0.1, 1, 4, 10)) +
  annotate("text", x = 1.8, y = 0.6, angle = 90, label = "Mean-time-to-capture\n1.67 years", size = 9, vjust = 1) +
  scale_y_continuous(labels = percent_format(), breaks = seq(0, 1, 0.05)) +
  few_mod() +
  ggtitle("Survival curve of running from the law") +
  ylab("Percent not caught") +
  xlab("Years from being put on the list") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 18),
        panel.grid.major.y = element_line(colour = "black", size = 0.1)) +
 #  geom_vline(xintercept = 613-92, colour = "red", size = 1, linetype = "dotted") +
  geom_vline(xintercept = 613 / 365.25, colour = "red", size = 1)
 #  geom_vline(xintercept = 613+92, colour = "red", size = 1, linetype = "dotted") 

# Take average time-to-capture of those captured
captured <- df[df$censor == 1, ]
mean(captured$days_to_capture, na.rm = T)

library(survival)
fit <- createSurvivalFrame(fit)

ggplot(fit[fit$time < 365 & fit$time >= 0, ]) + 
  geom_step(aes(x=time, y = surv, colour = strata), direction="hv") +
  geom_point(aes(x = time, y = surv, colour = strata))

ggplot(fi, aes(x = Freq)) +
  geom_density() +
  geom_vline(xintercept = quantile(tmp[tmp$Var1 == 1, "Freq"], 0.05), colour = "red") +
  geom_vline(xintercept = quantile(tmp[tmp$Var1 == 1, "Freq"], 0.95), colour = "red") +
  scale_x_continuous(breaks = seq(0, 1, 0.01), labels = percent_format())

library(muhaz)
fit2 <- muhaz(df$days_to_capture,
              df$censor,
              bw.method = "knn")
plot(fit2, xlim = c(0, 365.25),
     ylab = "Hazard Rate",
     xlab = "Days from being added to the list")
plot(fit2$pin$min.grid, fit2$bw.loc)
lines(fit2$est.grid, fit2$bw.loc.sm)
plot(fit2$bw.grid, fit2$globlmse)

plot(fit, fun = "log")

df$result_num <- as.integer(factor(df$result))

library(cmprsk)
ci <- cuminc(df$days_to_capture,
             df$censor,
             cencode = 1)
plot(ci, xlim = c(0, 365.25))


## Fit ######
library(flexsurv)
# install.packages("fitdistrplus")
library(fitdistrplus)
df$left <- df$days_to_capture
df$right <- ifelse(df$censor == 1, df$left, NA)

fit <- fitdistcens(censdata = df[df$left > 0, c("left", "right")], 
                   "weibull")
fit2 <- fitdistcens(censdata = df[df$left > 0, c("left", "right")], 
                   "lnorm")
fit3 <- fitdistcens(censdata = df[df$left > 0, c("left", "right")], 
                    "logis")
fit4 <- fitdistcens(censdata = df[df$left > 0, c("left", "right")], 
                    "norm")

cdfcompcens(list(fit, fit2, fit3, fit4),
            xlogscale = F,
            xlab = "Days since being added to the list\nnote: log10 scale",
            legendtext = c('Weibull fit', 'Log-normal fit', 'Logistic fit', 'Normal fit'),
            ylab = "Percentage caught",
            lines01 = T,
            lwd = 2.5)

df3 <- as.data.frame(list(x = seq(0.01, 0.99, 0.01),
          lognormal = unlist(quantile(fit2, seq(0.01, 0.99, 0.01))),
          raw_data = as.numeric(quantile(df$days_to_capture[df$censor == 1 & df$days_to_capture > 0], seq(0.01, 0.99, 0.01)))))

df3 <- melt(df3, "x")
ggplot(df3, aes(x = x, y = value, colour = factor(variable))) +
  geom_line(stat = "identity")
   


plotdistcens(df[df$left > 0, c("left", "right")],
             "lnorm",
             para = list(meanlog = fit2$estimate[1],
                         sdlog = fit2$estimate[2]))

df2 <- as.data.frame(list(times = qweibull(seq(0, 0.99, 0.1),
              shape = fit$estimate[1], scale = fit$estimate[2]),
     prob_of_being_caught = seq(0, 0.99, 0.1)))



ggplot(df2, aes(x = times / 10, y = prob_of_being_caught)) +
  geom_line() +
  geom_point() +
  coord_cartesian(xlim = c(0, 365.25))

fit <- fitdistcens(censdata = df[df$left >= 0 & df$right >= 0 & 
                                   !is.na(df$right) & !is.na(df$left), c("left", "right")], 
                   "lnorm")

plotdistcens(df[df$left >= 0 & df$right >= 0 & 
                  !is.na(df$right) & !is.na(df$left), c("left", "right")],
             "lnorm",
             para = list(meanlog = 19.8497178,
                         sdlog = 0.8706402))


fit <- fitdistcens(censdata = df[df$left >= 0 & df$right >= 0 & 
                                   !is.na(df$right) & !is.na(df$left), c("left", "right")], 
                   "gamma")

bfit <- bootdistcens(fit)

