---

#source('https://openmx.ssri.psu.edu/software/getOpenMx.R')
require(yarrr)
require(stargazer) #Hlavac, Marek (2018). stargazer: Well-Formatted Regression and Summary Statistics Tables.
# R package version 5.2.1. https://CRAN.R-project.org/package=stargazer
require(tidyverse)
require(OpenMx)
require(MASS)
require(psych) #for phi coeff
require(ggpubr)
require(kableExtra) #hmm, this does not work with Word output
require(flextable)
require(officer)
require(umx)
require(RVAideMemoire) #for bootstrapped CIs on Spearman correlations
require(ggExtra) #for scatterplot with marginals


source('http://ibg.colorado.edu/cdrom2016/verhulst/Power/powerFun.R')
#see http://www.people.vcu.edu/~bverhulst/power/powerScript.R
#mxOption(NULL, 'Default optimizer', 'NPSOL') #?? this gets 'could not find function mxOption'

tabnum <- 0 #keeps track of table numbers.
readprocessed <- 1 #if set to one, read processed data from OSF.
# If zero, recreates processed data file. 
# NB This option (0) will not work for external users, as it requires access to files with confidential data.
#---------------------------------------------------------------

if(readprocessed==0) {
  #---------------------------------------------------------------
  #Read data before any text as need to refer to Ns etc computed here.
  #Version from April 2019 includes LI indices using means rather than peaks, and also has measure of
  #cross-trial consistency
  #---------------------------------------------------------------
  mydir <-
    "/Users/dorothybishop/Dropbox/ERCadvanced/project twin kids/Project_files/"
  #mydir<-"c:/Users/pthompson/Dropbox/project twin kids/Project_files/"
  myfile <- "TwinsData_DATA_2019-04-12_1525.csv"
  mydata <- data.frame(read_csv(paste0(mydir, "Data/", myfile)))
  #str(mydata) #to inspect mydata
  w <-
    which(colnames(mydata) == 'include')  #this relates to exclusion criteria re IQ,autism etc : not used here
  colnames(mydata)[w] <-
    'myinclude' #to avoid problems with having 'include' as column name
}

if(readprocessed==0) {
  mydata$lat_index <-
    as.numeric(mydata$myli_mean)  #the original laterality_index values in file are LIs from peak
  #my_li_mean was added in April 2019
  #NB Reviewers want to see both analyses, so both retained
  
  mydata$myse_mean <- as.numeric(mydata$myse_mean)
  mydata$my_li_mean_even <- as.numeric(mydata$my_li_mean_even)
  mydata$my_li_mean_odd <- as.numeric(mydata$my_li_mean_odd)
  #Relevant column names. first block are from peak measures
  #print(colnames(mydata)[c(84:92, 101:105)]) - checking we have correct columns
  #rename columns to clarify which are mean and which are peak
  w <-
    which(
      colnames(mydata) %in% c(
        'laterality_index',
        'li_se',
        'low_ci',
        'high_ci',
        'lateralised_category'
      )
    )
  colnames(mydata)[w] <-
    c('peakLI',
      'peakLIse',
      'peaklow_ci',
      'peakhi_ci',
      'peak_latcat')
  
  #print(colnames(mydata)[c(84:92, 101:105)])
  
  mydata$qhp_handedness <- mydata$qhp_freq_r + .5 * mydata$qhp_freq_e
  
  #categorise laterality as 1=L, 2=bi, 3 = R depending on whether mean value is outside 1.96*SE limit
  #NB can use 1.96 or 1.65 here, but 1.96 is same as used by Wilson and Bishop
  CIterm <- 1.96
  mydata$li_cat <- NA
  w <- which(mydata$lat_index > 0)
  mydata$li_cat[w] <-
    1 #default is L laterlised for those with +ve values
  x <- mydata$lat_index[w] - CIterm * mydata$myse_mean[w]
  #need to be careful here re indices - we are now working with the w subset and finding which are -ve
  y <- which(x < 0)
  w1 <-
    w[y] #these should be the indices from the original full file with L sided lat but large SE -> bilateral
  mydata$li_cat[w1] <- 2
  
  #Now do the same for R lateralised
  w <- which(mydata$lat_index < 0)
  mydata$li_cat[w] <-
    3 #default is R lateralised for those with -ve values
  x <-
    mydata$lat_index[w] + CIterm * mydata$myse_mean[w] #this time we add the 1.96*SE term
  #need to be careful here re indices - we are now working with the w subset and finding which are -ve
  y <- which(x > 0)
  w1 <-
    w[y] #these should be the indices from the original full file with R sided lat but large SE -> bilateral
  mydata$li_cat[w1] <- 2
  
  #temp<-dplyr::select(mydata,lat_index,myse_mean,li_cat) #view temp to check output
  temp <-
    filter(mydata, n_trials > 11) #same as Wilson and bishop, have at least 12 trials
  oldnewtab <- table(temp$li_cat, temp$peak_latcat)
  #Mean method gives more bilateral cases - I guess not surprising.
  
  #--------------------------------------------------------------------------------------
  
}

if(readprocessed==0){
  #----------------------------------------------------------
  # Produce data frame of selected columns, called data.short
  #----------------------------------------------------------
  
  data.short<-dplyr::select(mydata,record_id,fam_id,age_at_test,female,zygosity,n_language_low,myinclude,twin,n_trials, lat_index,qhp_handedness,ehp_handedness,lang_probs,li_cat,ehp_handedness,qhp_handedness,ehp_write,peakLI,peakLIse,peak_latcat,mean_left_flow,mean_right_flow,mean_words,myse_mean,my_li_mean_even,my_li_mean_odd)
  
  #identify cases to exclude because of less than 12 trials on doppler
  w <- c(which(data.short$n_trials < 12), which(is.na(data.short$n_trials)))
  data.short$doppexcl <- 0
  data.short$doppexcl[w] <- 1
  w <- which(data.short$ehp_write == 999)
  data.short$ehp_write[w] <- NA
  data.short$lang_probs[data.short$lang_probs > 1] <- 1 #recodes so all values > 1 are now 1
  
  #-------------------------------------------------------
  # Create double entry file with twin 1 and 2 aligned
  #NB this code needs redoing to cope with possibility of additional columns
  # This works but depends crucially on the column numbers specified in the colnames column
  #-------------------------------------------------------
  nrec <- nrow(data.short)
  ncol <- ncol(data.short)
  Nexcluded <- vector()
  nuorder2 <- c(seq(from = 2, to = nrec, by = 2), seq(from = 1, to = nrec, by =
                                                        2))
  nuorder1 <- c(seq(from = 1, to = nrec, by = 2), seq(from = 2, to = nrec, by =
                                                        2))
  doubledata <- cbind(data.short[nuorder1, ], data.short[nuorder2, ])
  w1 <- which(colnames(doubledata) == 'record_id')
  
  colnames(doubledata)[w1[2]:ncol(doubledata)] <-
    paste0(colnames(doubledata)[1:(w1[2] - 1)], 2)
  #check all aligned
  check <- sum(doubledata$fam_id - doubledata$fam_id2)
  if (check > 0) {
    print('Twins not aligned!!!')
  }
  doubledatax <- filter(doubledata, myinclude > 0, myinclude2 > 0)
  Nexcluded[1] <- nrec - nrow(doubledatax)
  doubledatax <- filter(doubledata, n_trials > 11, n_trials2 > 11)
  Nexcluded[2] <- nrec - nrow(doubledatax) - Nexcluded[1]
  doubledatax <- filter(doubledatax, abs(peakLI) < 10, abs(peakLI2) < 10)
  doubledatax <-
    filter(doubledatax, abs(lat_index) < 10, abs(lat_index2) < 10)
  Nexcluded[3] <- nrec - nrow(doubledatax) - Nexcluded[1] - Nexcluded[2]
  
  #unit of analysis is twin pair excluded, so need to divide by 2 as one row per twin
  Nexcluded <- Nexcluded / 2
  
  #remove unwanted columns
  doubledata <- dplyr::select(doubledata,record_id, fam_id, age_at_test, female, zygosity, n_language_low , myinclude, twin, n_trials, lat_index, qhp_handedness , ehp_handedness , lang_probs, li_cat, ehp_write, peakLI,peak_latcat,myse_mean,doppexcl, mean_left_flow,mean_right_flow,mean_words,my_li_mean_even,my_li_mean_odd,female2, n_language_low2, myinclude2, twin2, n_trials2, lat_index2, qhp_handedness2, ehp_handedness2, lang_probs2, li_cat2, ehp_write2, peakLI2,peak_latcat2,mean_left_flow2,mean_right_flow2,mean_words2,my_li_mean_even,my_li_mean_odd,myse_mean2,doppexcl2)
  #make column where all DZ are 2 and all MZ are 1
  doubledata$MZDZ <- 1
  doubledata$MZDZ[doubledata$zygosity > 1] <- 2
  #make column for excluding pair if one or both has no good doppler data
  doubledata$doppout <- 0
  doubledata$doppout[doubledata$doppexcl == 1] <- 1
  doubledata$doppout[doubledata$doppexcl2 == 1] <- 1
  
  # Make langconcord column: this is just count of N twins with zero on lang_problems, else 1
  # This is not used for analysis but is reported in text
  doubledata$langconcord <- 0
  doubledata$langconcord[doubledata$lang_probs > 0] <- 1
  w <- which(doubledata$lang_probs2 > 0)
  doubledata$langconcord[w] <- 1 + doubledata$langconcord[w]
  
}

#-------------------------------------------------------------------------
if(readprocessed==1) {
  # Start here if reading in data from OSF
  doubledata <- read.csv('https://osf.io/m9jph/download') #when data public can read in direct from OSF
  #doubledata<-read.csv('TwinLatOSF.csv',stringsAsFactors = FALSE) #if reading downloaded file
  doubledata$fam_id <- doubledata$fam_randid #need fam_id variable at later stage
}
#-------------------------------------------------------------------------


#The twin concordance for language problems is reported as background information in the paper
conc.table <- table(doubledata$MZDZ, doubledata$langconcord, doubledata$twin)
conc.table <- conc.table[, , 1] #one twin pair only
rownames(conc.table) <- c('MZ', 'DZ')
colnames(conc.table) <- c('neither', 'one', 'both')
percenttab <- round(100 * prop.table(conc.table, 1), 0)
par(mfrow = c(1, 2))
par(mar = c(0, 0, 1, 0))
slices <- conc.table[1, 2:3]
lbls <- c('', '')
pie(slices, labels = lbls, main = "Language concordance: \nIdentical twins")
slices <- conc.table[2, 2:3]

pie(slices, labels = lbls, main = "Language concordance: \nNonidentical twins")
# Can view pie if wanted, but it's not used in the paper

#-------------------------------------------------------
# Exclude pairs with unusable data - either twin
#-------------------------------------------------------
doubledata$doppout <- doubledata$doppexcl + doubledata$doppexcl2
w <- which(doubledata$doppout == 2)
doubledata$doppout[w] <- 1

#remove case of extreme LI, defined as +/-10
w <- c(which(abs(doubledata$lat_index) > 12), which(abs(doubledata$lat_index2) >
                                                      12))
doubledata$doppout[w] <-1 #picks up case with LI of -21,which is more than 5 SD from mean!
# There seems no good reason for this extreme score- data look clean
# Nothing in this child's history to suggest abnormality. But LI is so extreme it is best removed.

doubledata2 <- filter(doubledata, doppout < 1)

# Unit of analysis is twin pair excluded, so need to divide by 2 as one row per twin

Npairwithhand <- nrow(doubledata) / 2
Npairwithdopp <- nrow(doubledata2) / 2

#identify N dropped because too few trials
w <- which(doubledata$n_trials < 12)
nlowN <- length(unique(doubledata$fam_id[w])) #This number will be reported in text


#Compute correlations which will be reported in the text
#check correlation between old (peak) and new (mean) LI methods
#use doubledata2 which just has included cases for Doppler
oldnewLIcor <-
  cor(doubledata2$peakLI, doubledata2$lat_index, use = 'pairwise.complete.obs')
oldnewLIcorDF <- nrow(doubledata2)- 2

#check split half reliablity
splithalfr <-
  cor(doubledata2$my_li_mean_even, doubledata2$my_li_mean_odd, use = 'pairwise.complete.obs')


#--------------------------------------------------------------
# Create table N twin pairs by zygo and sex
#--------------------------------------------------------------
tabnum<-tabnum+1
doubledata$zygosex<-10*doubledata$zygosity+doubledata$female
w<-which(doubledata$zygosex==31)
doubledata$zygosex[w] <-30
doubledata$zygosex<-as.factor(doubledata$zygosex)
levels(doubledata$zygosex) <- c('MZ female','MZ male','DZ female','DZ male','DZ male/female')
tab1 <- table(doubledata$zygosex,doubledata$doppout)
tab1f<-data.frame(cbind(tab1[1:5],tab1[1:5],tab1[6:10]))
tab1f[,1]<-tab1f[,2]+tab1f[,3]
tab1f[,1:2]<-tab1f[,1:2]/2
tab1f<-tab1f[,1:2]
tab1f<-data.frame(tab1f)
colnames(tab1f)<-c('All','With_fTCD')
nr<-nrow(tab1f)
tab1f[(nr+1),]<-colSums(tab1f)

forflex<-c(levels(doubledata$zygosex),'Total')
#stargazer(tab1f,type='text',title='Table 2: N twin pairs by zygosity and sex, whole sample and subset with fTCD',summary=FALSE)
mycaption<-paste0('Table ',tabnum,': N twin pairs by zygosity and sex')
#tab1f %>% kable(booktabs = T, caption = mycaption) %>% 
#  kable_styling()
tab1f<-cbind(forflex,tab1f)
colnames(tab1f)[1]<-'Twin type'
tab1f[,2]<-as.integer(tab1f[,2])
tab1f[,3]<-as.integer(tab1f[,3])
#calling flextable for first time gives huge N of Fontconfig warnings, but flextable
#does run. This seems a known problem ftp://ftp.oregonstate.edu/.1/cran/web/checks/check_results_flextable.html
ft1 <- flextable(tab1f)
ft1 <- width(ft1, width = 1.75)
ft1

#--------------------------------------------------------------


#Script for power analysis from Tim Bates
smallN <- nrow(doubledata2) / 2 #N pairs for laterality index; n= 141
bigN <- nrow(doubledata) / 2 #N pair for handedeness; n = 194

tmp = umx_make_TwinData(
  nMZpairs = 2000,
  nDZpairs = 2000,
  AA = .25,
  CC = 0,
  varNames = "var",
  mean = 0,
  empirical = TRUE
)
mzData = subset(tmp, zygosity == "MZ") #did not work, but substituted 2 lines below
dzData = subset(tmp, zygosity == "DZ")


ace = umxACE(
  selDVs = "var",
  sep = "_T",
  mzData = mzData,
  dzData = dzData,
  tryHard = "yes"
)
ae = umxModify(ace, "c_r1c1", name = "ae", tryHard = "yes")
nullModel = umxModify(ae, "a_r1c1", name = "e", tryHard = "yes")

#https://rdrr.io/cran/umx/man/power.ACE.test.html

#I had to get power by trial and error with the command below.
#We get to 80% power for aa of.215 with n = 139.

#The n is the total N twin pairs, MZ+DZ
thissize = .215 #use .185 for larger sample - with n = 189 this gives 80% power
p <-
  power.ACE.test(
    AA = thissize,
    CC = 0,
    update = "a_after_dropping_c",
    power = .73,
    method = "ncp"
  )

print(p)

thissize = .333 #
p <-
  power.ACE.test(
    AA = thissize,
    CC = 0,
    update = "a_after_dropping_c",
    power = .95,
    method = "ncp"
  )

print(p)




# function to do all density plots
makehistos <- function(mytitle){
  h<-ggdensity(doubledata, x = "thisx",
               add = "mean", rug = FALSE, title=mytitle,xlab='Score',ylab='Density',
               color = "Zygosity", fill = "Zygosity",palette = get_palette("aaas", 2)
  )
  return(h)
}

#function to compute Spearman correlation value
makespearman <- function(thismz,thisdz){
  #create text with Spearman correlation values.
  #NB need to compute these from *unjittered* data 
  MZrs <- cor(thismz$thisx,thismz$thisy, use="pairwise.complete.obs",method='spearman')
  MZn<-nrow(filter(thismz,!is.na(thisx),!is.na(thisy)))/2 #divide by 2 for N pairs
  DZrs <- cor(thisdz$thisx,thisdz$thisy, use="pairwise.complete.obs",method='spearman')
  DZn<-nrow(filter(thisdz,!is.na(thisx),!is.na(thisy)))/2 #divide by 2 for N pairs
  mylabelbit<-c(paste0('MZ: r (',(MZn-2),') = ',round(MZrs,3)),paste0('DZ: r (',(DZn-2),') = ',round(DZrs,3)))
  for (i in 1:2){
    while(nchar(mylabelbit[i])<18){
      mylabelbit[i]<-paste0(mylabelbit[i],'0')
    }
  }
  return(mylabelbit)
}

# function to do scatterplots
makescatter <- function(mytitle,xpos,ypos,filetitle,vline,mylabeldf,mylabelbit){
  g <-ggscatter(doubledata[1:singlerows,], x = 'myx', y = 'myy',
                xlab='Twin 1',ylab='Twin 2',
                shape = "Zygosity",
                # Change point shape by groups 
                color = "Zygosity", palette = get_palette("aaas", 2),
                # Color by groups 
                title = mytitle  )+
    geom_vline(xintercept=vline,linetype="dashed")+
    geom_hline(yintercept = vline,linetype="dashed") +
    annotate("text",label=mylabelbit[1], x = xpos[1],y = ypos[1],size=3)+
    annotate("text",label=mylabelbit[2], x = xpos[2],y = ypos[2],size=3)+
    
    
    ggsave(filetitle,width=3,height=3)
  return(g)
}



#------------------------------------------------------------------------
#scatterplot and correlation by zygosity as MZ/DZ: do for LI and handedness measures
#------------------------------------------------------------------------
#Use the makescatter function defined above for all the plots



doubledata$Zygosity <-
  as.factor(doubledata$MZDZ) #Initial caps for version that is factor with 2 values
levels(doubledata$Zygosity) <- c('MZ', 'DZ')
singlerows <- nrow(doubledata) / 2

#EHI plot
#compute spearman correlation using function and assemble into string
doubledata$thisx <- doubledata$ehp_handedness
doubledata$thisy <- doubledata$ehp_handedness2
thismz <- filter(doubledata, MZDZ == 1)
thisdz <- filter(doubledata, MZDZ == 2)
mylabelbit <- makespearman(thismz, thisdz)

#add jitter to points before plotting
myjitter <- rnorm(nrow(doubledata)) / 4.5 #need to add jitter to points for plotting
doubledata$myx <- doubledata$ehp_handedness + myjitter
myjitter <- rnorm(nrow(doubledata)) / 4.5
doubledata$myy <- doubledata$ehp_handedness2 + myjitter
mytitle <- 'Edinburgh Handedness Inventory'
filetitle <- 'ehiplot.png'
vline = 5.5
xpos <- c(2.2, 2.2)
ypos <- c(4.5, 3.5)
g1 <- makescatter(mytitle, xpos, ypos, filetitle, vline, mylabeldf, mylabelbit)
h1 <- makehistos(mytitle)


#compute spearman correlation using function and assemble into string
doubledata$thisx <- doubledata$qhp_handedness
doubledata$thisy <- doubledata$qhp_handedness2
thismz <- filter(doubledata, MZDZ == 1)
thisdz <- filter(doubledata, MZDZ == 2)
mylabelbit <- makespearman(thismz, thisdz)

myjitter <- rnorm(nrow(doubledata)) / 4
doubledata$myx <- doubledata$qhp_handedness + myjitter
myjitter <- rnorm(nrow(doubledata)) / 4
doubledata$myy <- doubledata$qhp_handedness2 + myjitter
xpos <- c(3, 3) #coordinates for correlation on plot
ypos <- c(8, 6)
mytitle <- 'Quantification of Hand Preference'
filetitle <- 'qhpplot.png'
vline = 11
g2 <- makescatter(mytitle, xpos, ypos, filetitle, vline, mylabeldf, mylabelbit)
h2 <- makehistos(mytitle)

w <- which(doubledata$doppout == 1) #rather than using doubledata2, we just create column which
#turns LIs to NA for those with too few trials
doubledata$lat_index[w] <- NA 
doubledata$lat_index2[w] <- NA
doubledata$peakLI[w]<-NA
doubledata$peakLI2[w] <- NA

#compute spearman correlation using function and assemble into string
doubledata$thisx <- doubledata$lat_index
doubledata$thisy <- doubledata$lat_index2
thismz <- filter(doubledata, MZDZ == 1)
thisdz <- filter(doubledata, MZDZ == 2)
mylabelbit <- makespearman(thismz, thisdz)
doubledata$myx <- doubledata$lat_index #no need for jitter with LI
doubledata$myy <- doubledata$lat_index2
mytitle <- 'Laterality Index (fTCD)'
xpos <- c(4, 4) #coordinates for correlation on plot
ypos <- c(10, 9)
filetitle <- 'LIplot.png'
vline = 0
g3 <-
  makescatter(mytitle, xpos, ypos, filetitle, vline, mylabeldf, mylabelbit)
h3 <- makehistos(mytitle)

#add histogram for peak measure, as requested by reviewers
doubledata$thisx <- doubledata$peakLI
doubledata$thisy <- doubledata$peakLI2
w <-
  which(doubledata$doppout == 1) #rather than using doubledata2, we just create column which
#turns LIs to NA for those with too few trials
doubledata$peakindex[w] <- NA #this removes outlier value
doubledata$peakindex2[w] <- NA
thismz <- filter(doubledata, MZDZ == 1)
thisdz <- filter(doubledata, MZDZ == 2)
mytitle <- 'Laterality Index (peak fTCD)'

filetitle<-'LIplotpeak.png'
vline=0
h4<-makehistos(mytitle)
ggsave("peakLI.png",width=4,height=4,dpi=300)#save most recent plot
###########################################################################################

# Make one big plot with each task in a panel. Save as EPS format.
ggarrange(g1,g2,g3,ncol=1,nrow=3,common.legend = TRUE,heights=c(1,1,1))
ggsave("combined_plot.png",width=4,height=12,dpi=300)

ggarrange(h1,h2,h3,ncol=1,nrow=3,common.legend = TRUE,heights=c(1,1,1))
ggsave("combined_densities.png",width=4,height=12,dpi=300)

setEPS()
postscript("combined_plot.eps",width=4,height=12)
ggarrange(g1,g2,g3,ncol=1,nrow=3,common.legend = TRUE,heights=c(1,1,1))
dev.off()

setEPS()
postscript("combined_densities.eps",width=4,height=12)
ggarrange(h1,h2,h3,ncol=1,nrow=3,common.legend = TRUE,heights=c(1,1,1))
dev.off()


dopoly<-1
if(dopoly==1){
  #AE model; try ordinal version for handedness.
  # Based on ibg.colorado.edu/cdrom2010/medland/wed_morning/wed_morning2010.pdf
  
  
  #use double-entry file ? otherwise problems if null value for one category
  myehp<-dplyr::select(doubledata,MZDZ,ehp_handedness,ehp_handedness2)
  npair<-nrow(myehp)
  my.mz<-filter(myehp[1:npair,],MZDZ==1)
  mypolytab.mz<-table(my.mz$ehp_handedness,my.mz$ehp_handedness2)
  my.dz<-filter(myehp[1:npair,],MZDZ==2)
  mypolytab.dz<-table(my.dz$ehp_handedness,my.dz$ehp_handedness2)
  Emypoly.mz<-polychoric(mypolytab.mz)$rho
  Emypoly.dz<-polychoric(mypolytab.dz)$rho
  print('Polychoric correlations MZ/DZ (rho)')
  print(c(Emypoly.mz,Emypoly.dz))
  
  #QHP has too many values for computing polychoric, so reduced to broader categories by dividing by 4.
  #This is arbitrary - choice of divisor does make a slight difference
  myqhp<-dplyr::select(doubledata,MZDZ,qhp_handedness,qhp_handedness2)
  npair<-nrow(myqhp)
  my.mz<-filter(myqhp[1:npair,],MZDZ==1)
  mypolytab.mz<-table(round(my.mz$qhp_handedness/4,0),round(my.mz$qhp_handedness2/4,0))
  my.dz<-filter(myqhp[1:npair,],MZDZ==2)
  mypolytab.dz<-table(round(my.dz$qhp_handedness/4,0),round(my.dz$qhp_handedness2/4,0))
  Qmypoly.mz<-polychoric(mypolytab.mz)$rho
  Qmypoly.dz<-polychoric(mypolytab.dz)$rho
  Qmypoly.mz
  Qmypoly.dz
  
}

#--------------------------------------------------------------
#It's clear handedness all measures are v non-normal.
#--------------------------------------------------------------
donorm<-0
if (donorm==1){
  shapiro.test(doubledata$ehp_handedness)
  qqnorm(doubledata$ehp_handedness)
  shapiro.test(doubledata$qhp_handedness)
  qqnorm(doubledata$qhp_handedness)
  shapiro.test(doubledata$lat_index)
  qqnorm(doubledata$lat_index)
  shapiro.test(doubledata$peakLI)
  qqnorm(doubledata$peakLI) #normality test added for peak LI
  shapiro.test(doubledata$taskssum)
  qqnorm(doubledata$taskssum)
  hist(doubledata$ehp_handedness,breaks=20)
  hist(doubledata$qhp_handedness,breaks=20)
  hist(doubledata$lat_index,breaks=20)
}

#A bit more data wrangling to create measures that are suggested by reviewer comments

#We'll compare with Somers et al, so define a binary category from the peak LI
doubledata$peakcatbinary <- doubledata$peak_latcat
w <- which(doubledata$peakcatbinary < 0)
doubledata$peakcatbinary[w] <- 0 #combines bilateral and R
doubledata$peakcatbinary2 <- doubledata$peak_latcat2
w <- which(doubledata$peakcatbinary2 < 0)
doubledata$peakcatbinary2[w] <- 0 #combines bilateral and R





#---------------------------------------------------------------
#Try fitting AE model, despite low correlations.
#--------------------------------------------------------------
# Predict it won't fit because of zero DZ correlation - only compatible with no genetic effect! 
# Surprisingly, fit seems OK.
## NB However, really need a model that takes into account nonnormality! - see below re ordinal data
#Retain this because it is needed to define aatable
nv <- 1 #N variables
ntv <- nv*2 #N columns (twins)
aatable <-data.frame(matrix(NA,nrow=7,ncol=6))
colnames(aatable)<-c('Measure','N_MZ','N_DZ','a','a_SE','a_sqr')
mxcols <-matrix(c('ehp_handedness','ehp_handedness2','qhp_handedness','qhp_handedness2','lat_index','lat_index2','peakLI','peakLI2','peak_latcat','peak_latcat2','mean_left_flow','mean_left_flow2','mean_right_flow','mean_right_flow2'),nrow=7,byrow=TRUE)
for (mytask in 1:5){
  selVars <- mxcols[mytask,1:2]
  latcols<-which(names(doubledata)%in%selVars)
  rowrange<-1:(nrow(doubledata)/2)
  mymz <- filter(doubledata[rowrange,],MZDZ==1)
  mydz <- filter(doubledata[rowrange,],MZDZ==2)
  if(mytask>2)
  {mymz<-filter(mymz,doppout==0)
  mydz<-filter(mydz,doppout==0)
  }
  MZdata  <-  mxData(scale(mymz[,latcols]), type="raw" )
  DZdata  <-  mxData(scale(mydz[,latcols]), type="raw" )
  a <-  mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.6, label="a11", name="a" )
  e <-  mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=.6, label="e11", name="e" )
  A <-  mxAlgebra( expression=a %*% t(a), name="A" )
  E <-  mxAlgebra( expression=e %*% t(e), name="E" )
  Mean    <-  mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values= 0, label="mean", name="Mean" )
  expMean <-  mxAlgebra( expression= cbind(Mean,Mean), name="expMean")
  expCovMZ <- mxAlgebra( expression= rbind  ( cbind(A+E , A),cbind(A   , A+E)), name="expCovMZ" )
  expCovDZ <- mxAlgebra( expression= rbind  ( cbind(A+E     , 0.5%x%A),cbind(0.5%x%A , A+E)),  name="expCovDZ" ) 
  obs <-  list(a,e,A,E,Mean,expMean,expCovMZ,expCovDZ)
  fun <- mxFitFunctionML()
  mzExp <- mxExpectationNormal(covariance="expCovMZ", means="expMean", dimnames=selVars )
  dzExp <- mxExpectationNormal(covariance="expCovDZ", means="expMean", dimnames=selVars )
  MZ <- mxModel("MZ", obs, MZdata, fun, mzExp)
  DZ <- mxModel("DZ", obs, DZdata, fun, dzExp)
  aeFun <- mxFitFunctionMultigroup(c("MZ","DZ"))
  ae <- mxModel("AE", MZ, DZ, fun, aeFun)
  aeFit <- mxRun(ae, silent = T)
  summary(aeFit)
  pars <- round(cbind(aeFit$output$estimate, aeFit$output$standardErrors), 3)
  colnames(pars) <- c("Estimates", "Std. Err.")
  aatable[mytask,1]<-mxcols[mytask,1]
  aatable[mytask,2]<-nrow(mymz)
  aatable[mytask,3]<-nrow(mydz)
  aatable[mytask,4]<-pars[1,1]
  aatable[mytask,5]<-pars[1,2]
  aatable[mytask,6]<-round(pars[1,1]^2,2)
  # aatable has the a squared estimates that we need in final column, but the SE is for unsquared a!
  # We want a SE for heritability estimate, ie squared value. This will be achieved by bootstrapping.
  write.csv(aatable,'aatable.csv',sep='',row.names=F,col.names=T)
}

#Including formatting into a nice table with everything in it
# uutable includes bootstrapped confidence intervals

# NB these aren't computed for the flow measures, as the AE model is not suitable for these!
# Flow measures look as if they have influence of CE, so this is computed for them!

mycols<-c("ehp_handedness", "qhp_handedness", "lat_index","peakLI","peakcatbinary","mean_left_flow","mean_right_flow")
mycolnames<-c("Edinburgh Handedness Inventory","Quantification of Hand Preference","Laterality index (fTCD)","Laterality index (peak FTCD)", "Binary typical (1)/atypical (0)","Mean flow L","Mean flow R")

uutable<-data.frame(matrix(NA,nrow=(length(mycols)),ncol=10))
colnames(uutable)<-c('Measure','N MZ pair','N DZ pair','rMZ','rDZ','a2','chisq','p','RMSEA')
#NB 95CI will be done from bootstrap_CI.Rmd

umx_set_auto_plot(FALSE) #turn off the umx plots
#umx gives v slightly different estimates of a2 than the standard ACE above.
#I am using the OpenMx versions 


df = doubledata # for ease of typing
#create new columns, X_T1 and X_T2, which are recycled for each measure
for (j in 1:length(mycols)){
  
  if (j>2)
    # {df<-filter(df,myse_mean<1,doppout==0)}
  {df<-subset(df,doppout==0)} #drop rows with missing data for the doppler measures
  
  col1<-which(names(df)==mycols[j]) #find the column for this measure t1
  df$X_T1 <-df[,col1] #allocate contents of this column to variable X_T1
  name2<-paste0(mycols[j],'2') #repeat for twin 2
  col2<-which(names(df)==name2)
  df$X_T2 <-df[,col2]
  
  mzData = subset(df, MZDZ==1)
  dzData = subset(df, MZDZ==2)
  
  #Populate first cols of uutable
  thisrow<-j
  uutable$Measure[thisrow] <- mycolnames[j]
  uutable$`N MZ pair`[thisrow] <- nrow(mzData)/2
  uutable$`N DZ pair`[thisrow] <- nrow(dzData)/2
  sp.m<-spearman.ci(mzData$X_T1, mzData$X_T2, nrep = 1000, conf.level = 0.95)
  uutable$rMZ[thisrow] <- paste0(round(sp.m$estimate,2)," (",round(sp.m$conf.int[1],2),' to ',round(sp.m$conf.int[2],2),')')
  sp.d<-spearman.ci(dzData$X_T1, dzData$X_T2, nrep = 1000, conf.level = 0.95)
  uutable$rDZ[thisrow] <- paste0(round(sp.d$estimate,2)," (",round(sp.d$conf.int[1],2),' to ',round(sp.d$conf.int[2],2),')')
  
  # Test for genetic influence 
  m1 = umxACE("EHP", selDVs = "X", sep="_T", mzData = mzData, dzData= dzData)
  m2 = umxModify(m1, "c_r1c1", name="ae", comparison= TRUE)
  m3 = umxModify(m2, "a_r1c1", name="e", comparison= TRUE)
  m4 = umxModify(m1, "a_r1c1", name="ce", comparison= TRUE)
  mycomp<-umxCompare(m2, m3)
  summ<-summary(m2)
  if(j>5){
    mycomp<-umxCompare(m4, m3) #CE model for the flow measures
    summ<-summary(m4)
  }
  
  #compute a2 as proportion of variance from raw estimates
  a2<-summ$parameters$Estimate[2]^2/(summ$parameters$Estimate[2]^2+summ$parameters$Estimate[3]^2)
  if (sp.d$estimate>sp.m$estimate){a2<-a2*-1}
  uutable$a2[thisrow]<-round(a2,2)
  
  #add chi sq for comparison with model 3
  uutable$chisq[thisrow]<-round(mycomp$`∆ -2LL`[2],1)
  uutable$p[thisrow]<-mycomp$p[2]
  uutable$RMSEA[thisrow]<-summ$RMSEA
  if(j<6){
    #add bootstrapped CI computed with bootstrap_ci.Rmd
    
    #########################################################################
    #need to specify directory to find bootstrapped estimates in filepath
    #########################################################################
    
    # bootfile<-paste0('data/',mycols[j],'_bootstraps_umx_5000.csv')
    # myboots<-read.csv(bootfile,stringsAsFactors=F)
    # goodresults<-myboots[myboots$infoDefinite==T,]
    # q3<-quantile(goodresults$apropall,.05)
    # q4<-quantile(goodresults$apropall,.95)
    # myCI<-(paste(round(q3,2),' to ',round(q4,2)))
    # uutable$a2[j]<-paste0(uutable$a2[j],"(",myCI,")")
  }
}

#write.csv(uutable,'uutable.csv',sep='',row.names=F,col.names=T)




#Multivariate with both handedness measures
#Not included in paper: but requested by reviewer and added in response
multivar<-1
if(multivar==1){
  umx_set_auto_plot(TRUE)
  #reformat data to fit umx format
  mytemp<-dplyr::select(doubledata2,ehp_handedness,qhp_handedness,lat_index,ehp_handedness2,qhp_handedness2,lat_index2,MZDZ)
  colnames(mytemp)[1:6]<-c('EHP_1','QHP_1','Lat_index_1','EHP_2','QHP_2','Lat_index_2')
  mzData<-filter(mytemp,MZDZ==1)
  dzData<-filter(mytemp,MZDZ==2)
  
  m1 = umxACE(selDVs = c("EHP", "QHP","Lat_index"), sep="_", mzData = mzData, dzData= dzData)
  summary(m1) #confirms that all C paths are zero
  
  m2 = umxModify(m1, regex = "c_r.c[1-3]", name="dropC_123", comparison= TRUE)
  umxSummary(m2) #displays path diagram
  summary(m2)
  mxCompare(m1,m2) #gives p value of 1!  Dropping C for all 3 vars has no effect
  write.csv(summary(m2)$parameters,'summarym2.csv')
  
  #now drop A terms 
  #  m3 = umxModify(m2, regex = "a_r.c[1-3]", name="dropA_123", comparison= TRUE)
  #  umxSummary(m3) #displays path diagram 
  # umxCompare(m2, m3)
  # mxCompare(m2,m3) #p = .107, i.e. can drop A terms as well without sig loss of fit!
  # }
  
  #Extract parameters from model 2 for report. These need to be converted to standardized format.  Can check against the paths shown in the figure.
  
  # First find the sum of estimates for all paths leading to each measure.
  
  myparam <- summary(m2)$parameters
  myparam<-myparam[-(c(1:3)),]
  totsq<-c(0,0,0) #total sum squared estimates for each of the 3 measures
  totsq[1] <- sum(myparam$Estimate[myparam$row==1]^2)
  totsq[2] <- sum(myparam$Estimate[myparam$row==2]^2)
  totsq[3] <- sum(myparam$Estimate[myparam$row==3]^2)
  myparam$row<-as.numeric(myparam$row)
  myparam$col<-as.numeric(myparam$col)
  myparam$standardised<-NA
  for (i in 1:nrow(myparam)){
    myparam$standardised[i] <- sqrt(myparam$Estimate[i]^2/totsq[myparam$row[i]])
  }
  myparam$s.se<-myparam$Std.Error*myparam$standardised/myparam$Estimate
  
  #ensure sign is correct
  w<-which(myparam$Estimate<0)
  myparam$standardised[w]<--1*myparam$standardised[w]
  
  #Now put in a table
  multdf <- data.frame(matrix(ncol = 8, nrow = 6))
  names(multdf)<-c("Term", "Measure", "Estimate.EQL","SE.EQL","Estimate.QL","SE.QL","Estimate.L","SE.L")
  multdf$Term <-c('A','A','A','E','E','E')
  multdf$Measure <- c('EHI','QHP','LI','EHI','QHP','LI')
  multdf[,3]<-round(myparam$standardised[myparam$col==1],3)
  multdf[,4]<-abs(round(myparam$s.se[myparam$col==1],3))
  multdf[c(2,3,5,6),5]<-round(myparam$standardised[myparam$col==2],3)
  multdf[c(2,3,5,6),6]<-abs(round(myparam$s.se[myparam$col==2],3))
  multdf[c(3,6),7]<-round(myparam$standardised[myparam$col==3],3)
  multdf[c(3,6),8]<-abs(round(myparam$s.se[myparam$col==3],3))
  
  write.csv(multdf,'multdftab.csv')
}


#stargazer(aatable,type='text',title='Table 3: Heritability estimates from AE model for different measures: ignoring nonnormality',summary=FALSE)
tabnum<-tabnum+1

mycaption<-paste0('Table ',tabnum,': Heritability estimates from AE model for different measures: ignoring nonnormality')
#aatable %>% kable(booktabs = T, caption = mycaption) %>% 
#  kable_styling()

ft <-flextable(data = uutable[1:5,1:9])

ft <-width(ft, j=1,width = 2.5) #sets width of column j to 2.5 inches
ft

#NB all tasks have been scaled; skew in the handedness measures has not been accounted for

doumx <- function(col1,col2){
  mynrow<-nrow(doubledata)
  zygo<-doubledata$MZDZ
  myhand1<-round(doubledata[,col1],0) #just to fit labelling convention
  myhand2<-round(doubledata[,col2],0)
  for.umx<-as.data.frame(cbind(zygo,myhand1,myhand2))
  colnames(for.umx)<-c('zygo','hand1','hand2')
  for.umx$hand1<-as.factor(for.umx$hand1)
  for.umx$hand2<-as.factor(for.umx$hand2)
  nlevels<-nlevels(for.umx$hand1)
  
  for.umx<-for.umx[1:(mynrow/2),] #make single entry
  
  ordDVs = c("hand1", "hand2")
  selDVs = c("hand")
  # Make the ordinal variables into mxFactors (ensure ordered is TRUE, and require levels)
  
  for.umx[, ordDVs] <- mxFactor(for.umx[, ordDVs], levels = c(0:nlevels))
  
  mzData <- for.umx[for.umx$zygo ==1, ]
  dzData <- for.umx[for.umx$zygo ==2, ]
  m1 = umxACE(selDVs = selDVs, dzData = dzData, mzData = mzData, sep = '')
  m2 = umxModify(m1, update="c_r1c1",name = "drop_C")
  m3 = umxModify(m2, "a_r1c1", name="drop_E", comparison= TRUE)
  umxCompare(m2, m3)
  umxSummary(m2)
  
  umxSummary(m2,showEstimates='std',SE=TRUE,report='markdown',RMSEA_CI=TRUE)
  aSE<-m2$output$standardErrors[nlevels] #nlevels does identify the correct row here
  Qa<-m2$output$estimate[nlevels]
  QlowCI<-(Qa-1.96*aSE)
  QhiCI<-(Qa+1.96*aSE)
  umxresult<-as.data.frame(c(Qa,aSE,QlowCI,QhiCI,Qa^2))
  
  
  return(umxresult)
}

# ===================
# = Ordinal example from Bates with handedness data substituted=
# ===================
umxresult<-data.frame(matrix(NA,nrow=5,ncol=2))
colnames(umxresult)<-c('EHI','QHP')
rownames(umxresult)<-c('a','a.SE','aLowCI','aHighCI','a2')
for (i in 1:2){
  mycolname<-'ehp_handedness'
  if (i==2){
    mycolname<-'qhp_handedness'
  }
  col1<-which(colnames(doubledata)==mycolname)
  col2<-which(colnames(doubledata)==paste0(mycolname,'2'))
  
  umxresult[,i]<-doumx(col1,col2)
  
}


dftable <- data.frame(matrix(nrow = 3, ncol = 10))
colnames(dftable)<-c('Measure','N MZ','N DZ','MZ proband mean (SD)','DZ proband mean (SD)',
                     'MZ co-twin mean (SD)','DZ co-twin mean (SD)','t','df','p')
dftable[,1]<-aatable[3,1] #names of measures

# T-tests as for DeFries Fulker analysis
Ecutoff<-6
Elefties<-filter(doubledata,ehp_handedness<Ecutoff)
dftable[1,2]<-length(which(Elefties$MZDZ==1))
dftable[1,3]<-length(which(Elefties$MZDZ==2))
t.test(Elefties$ehp_handedness~Elefties$MZDZ)
myt<-t.test(Elefties$ehp_handedness2~Elefties$MZDZ)
aggmeans <-aggregate(Elefties$ehp_handedness, by=list(Elefties$MZDZ),
                     FUN=mean, na.rm=TRUE)
aggsds <-aggregate(Elefties$ehp_handedness, by=list(Elefties$MZDZ),
                   FUN=sd, na.rm=TRUE)

dftable[1,4]<-paste0(round(aggmeans[1,2],2)," (",round(aggsds[1,2],2),")")
dftable[1,5]<-paste0(round(aggmeans[2,2],2)," (",round(aggsds[2,2],2),")")
aggmeans <-aggregate(Elefties$ehp_handedness2, by=list(Elefties$MZDZ),
                     FUN=mean, na.rm=TRUE)
aggsds <-aggregate(Elefties$ehp_handedness2, by=list(Elefties$MZDZ),
                   FUN=sd, na.rm=TRUE)
dftable[1,6]<-paste0(round(aggmeans[1,2],2)," (",round(aggsds[1,2],2),")")
dftable[1,7]<-paste0(round(aggmeans[2,2],2)," (",round(aggsds[2,2],2),")")
dftable[1,8]<-round(myt$statistic,2)
dftable[1,9]<-round(myt$parameter,1)
dftable[1,10]<-round(myt$p.value,3)

Qcutoff<-11
Qlefties<-filter(doubledata,qhp_handedness<Qcutoff)
thisrow<-2
dftable[thisrow,2]<-length(which(Qlefties$MZDZ==1))
dftable[thisrow,3]<-length(which(Qlefties$MZDZ==2))
t.test(Qlefties$qhp_handedness~Qlefties$MZDZ)
myt<-t.test(Qlefties$qhp_handedness2~Qlefties$MZDZ)
aggmeans <-aggregate(Qlefties$qhp_handedness, by=list(Qlefties$MZDZ),
                     FUN=mean, na.rm=TRUE)
aggsds <-aggregate(Qlefties$ehp_handedness, by=list(Qlefties$MZDZ),
                   FUN=sd, na.rm=TRUE)

dftable[thisrow,4]<-paste0(round(aggmeans[1,2],2)," (",round(aggsds[1,2],2),")")
dftable[thisrow,5]<-paste0(round(aggmeans[2,2],2)," (",round(aggsds[2,2],2),")")
aggmeans <-aggregate(Qlefties$qhp_handedness2, by=list(Qlefties$MZDZ),
                     FUN=mean, na.rm=TRUE)
aggsds <-aggregate(Qlefties$qhp_handedness2, by=list(Qlefties$MZDZ),
                   FUN=sd, na.rm=TRUE)
dftable[thisrow,6]<-paste0(round(aggmeans[1,2],2)," (",round(aggsds[1,2],2),")")
dftable[thisrow,7]<-paste0(round(aggmeans[2,2],2)," (",round(aggsds[2,2],2),")")
dftable[thisrow,8]<-round(myt$statistic,2)
dftable[thisrow,9]<-round(myt$parameter,1)
dftable[thisrow,10]<-round(myt$p.value,3)

LIcutoff<-0
Drighties<-filter(doubledata,lat_index<LIcutoff)
thisrow<-3
dftable[thisrow,2]<-length(which(Drighties$MZDZ==1))
dftable[thisrow,3]<-length(which(Drighties$MZDZ==2))
t.test(Drighties$lat_index~Drighties$MZDZ)
myt<-t.test(Drighties$lat_index2~Drighties$MZDZ)
aggmeans <-aggregate(Drighties$lat_index, by=list(Drighties$MZDZ),
                     FUN=mean, na.rm=TRUE)
aggsds <-aggregate(Drighties$ehp_handedness, by=list(Drighties$MZDZ),
                   FUN=sd, na.rm=TRUE)

dftable[thisrow,4]<-paste0(round(aggmeans[1,2],2)," (",round(aggsds[1,2],2),")")
dftable[thisrow,5]<-paste0(round(aggmeans[2,2],2)," (",round(aggsds[2,2],2),")")
aggmeans <-aggregate(Drighties$lat_index2, by=list(Drighties$MZDZ),
                     FUN=mean, na.rm=TRUE)
aggsds <-aggregate(Drighties$lat_index2, by=list(Drighties$MZDZ),
                   FUN=sd, na.rm=TRUE)
dftable[thisrow,6]<-paste0(round(aggmeans[1,2],2)," (",round(aggsds[1,2],2),")")
dftable[thisrow,7]<-paste0(round(aggmeans[2,2],2)," (",round(aggsds[2,2],2),")")
dftable[thisrow,8]<-round(myt$statistic,2)
dftable[thisrow,9]<-round(myt$parameter,1)
dftable[thisrow,10]<-round(myt$p.value,3)

#just for interest: try regression analysis with probands selected by R sided language
Dcat<-filter(doubledata,li_cat>2)
dfreg1<-lm(Dcat$lat_index2~Dcat$lat_index*Dcat$MZDZ)
summary(dfreg1)

tabnum<-tabnum+1
mycaption<-paste0('Table ',tabnum,': Extreme probands, t-tests comparing MZ/DZ co-twins')
#dftable %>% kable(booktabs = T, caption = mycaption) %>% 
#  kable_styling()
tdftable<-data.frame(t(dftable))
colnames(tdftable)<-c('Edinburgh Handedness','QHP','Language LI')
tdftable<-tdftable[-1,] #top row will now be header

rowvec<-rownames(tdftable)
tdftable<-cbind(rowvec,tdftable)
colnames(tdftable)[1]<-' '
ft3<-flextable(tdftable)
ft3<-width(ft3, width = 1.75)
ft3
#https://davidgohel.github.io/flextable/articles/examples.html

if(readprocessed==0){
  mylm <- summary(lm(doubledata$mean_left_flow ~ doubledata$age_at_test+doubledata$female+doubledata$mean_words))
  doubledata$Lflow.resid<-doubledata$mean_left_flow -(mylm$coefficients[1,1]+
                                                        mylm$coefficients[2,1]*doubledata$age_at_test+
                                                        mylm$coefficients[3,1]*doubledata$female)
  doubledata$Lflow.resid2<-doubledata$mean_left_flow2 -(mylm$coefficients[1,1]+
                                                          mylm$coefficients[2,1]*doubledata$age_at_test+
                                                          mylm$coefficients[3,1]*doubledata$female)
  
  mylm2 <- summary(lm(doubledata$mean_right_flow ~ doubledata$age_at_test+doubledata$female+doubledata$mean_words))
  doubledata$Rflow.resid<-doubledata$mean_right_flow -(mylm2$coefficients[1,1]+
                                                         mylm2$coefficients[2,1]*doubledata$age_at_test+
                                                         mylm2$coefficients[3,1]*doubledata$female)
  doubledata$Rflow.resid2<-doubledata$mean_right_flow2 -(mylm2$coefficients[1,1]+
                                                           mylm2$coefficients[2,1]*doubledata$age_at_test+
                                                           mylm2$coefficients[3,1]*doubledata$female)
}                                                   

#change the ID number - make random set of nums
if(readprocessed==0){
  npairs<-nrow(doubledata)/2
  nupairnum<-sample(1000, npairs, replace = FALSE)
  doubledata$rand_id[1:npairs]<-paste0(nupairnum,'A')
  doubledata$rand_id[(npairs+1):(2*npairs)]<-paste0(nupairnum,'B')
  doubledata$fam_randid<-rep(nupairnum,2)
  dataforOSF<-dplyr::select(doubledata,c(rand_id,fam_randid,age_at_test,female,zygosity,n_language_low,myinclude,twin,n_trials,lat_index,myse_mean,qhp_handedness,ehp_handedness,lang_probs,li_cat,ehp_write,peakLI,peak_latcat,mean_left_flow,mean_right_flow,Lflow.resid,Rflow.resid,mean_words,doppexcl,female2,n_language_low2,myinclude2,twin2,n_trials2,lat_index2,myse_mean2,qhp_handedness2,ehp_handedness2,lang_probs2,li_cat2,ehp_write2,peakLI2,peak_latcat2,mean_left_flow2,mean_right_flow2,Lflow.resid2,Rflow.resid2,doppexcl2,MZDZ,langconcord,doppout,zygosex,my_li_mean_even,my_li_mean_odd))
  dataforOSF$age_at_test<-round(dataforOSF$age_at_test/24,0) #age in 2 yr bands to improve anonymisation
  dataforOSF$age_at_test<-as.factor(dataforOSF$age_at_test)
  levels(dataforOSF$age_at_test)<-c('6_7','8_9','10_11','12_13')
  write.csv(dataforOSF,'TwinLatOSF.csv',row.names=FALSE)
}

#added analysis of handedness by LI, using peak_LI divided at zero and writing hand
#
doubledata2$peakbinaryLI <-0
doubledata2$peakbinaryLI[doubledata2$peakLI>0]<-1
doubledata2$binaryLI <-0
doubledata2$binaryLI[doubledata2$lat_index>0]<-1
peakhandLItab <-table(doubledata2$peakbinaryLI,doubledata2$ehp_write)
peakphandLItab <-prop.table(peakhandLItab, 2) # row percentages 

handLItab <-table(doubledata2$binaryLI,doubledata2$ehp_write)
phandLItab <-prop.table(handLItab, 2) # row percentages
#proportion test against Carey figures
#L handers
p0 <- phandLItab[2,1]
pE <- .67
q <- 1-p0
n <- handLItab[1,1]+handLItab[2,1]
zL<- (p0-pE)/sqrt(p0*q/n)

#R handers
p0 <- phandLItab[2,2]
pE <- .91
q <- 1-p0
n <- handLItab[2,1]+handLItab[2,2]
zR<- (p0-pE)/sqrt(p0*q/n)

#Requested by reviewer CM

# classic plot :
p <- ggplot(doubledata2, aes(x=lat_index, y=peakLI,color=as.factor(ehp_write))) +
geom_point(shape = 18 ) + 
labs(x = "Laterality index (mean)",y="Laterality index (peak)")+
labs(colour = "R handed")+
theme(legend.position="bottom")

# marginal density
p2 <- ggMarginal(p, type="density")
ggsave("peak_vs_meanLI.png",plot=p2,width=4,height=4,dpi=300)




