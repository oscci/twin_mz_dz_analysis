require(umx)
# Double checking power calculations, after query by Guy Vingerhoets.
# Using script from https://cran.r-project.org/web/packages/umx/index.html, see p 23.

# NOTE: the value of AA entered into this function is proportion of variance explained by
# genetic factors in AE model (after dropping c), i.e equivalent to standardized a^2.
# Our sample Ns are:
# twins1   = c(96, 98) # N MZ/DZ twin pairs, EHI
# twins2   = c(65, 76) # N twin pairs, LI

#power.ACE.test(AA= .2, CC= 0, update= "a_after_dropping_c", MZ_DZ_ratio= .97, n= 194, power= NULL, method= "empirical")
#power.ACE.test(AA= .3, CC= 0, update= "a_after_dropping_c", MZ_DZ_ratio= .97, n= 194, power= NULL, method= "empirical")

#power.ACE.test(AA= .2, CC= 0, update= "a_after_dropping_c", MZ_DZ_ratio= .85, n= 141, power= NULL, method= "empirical")


# ==========================================================================================
# = Make vectors of A and power values a loop over these to build an ncp-based table of Ns =
# ==========================================================================================

myAA     = seq(.15,.35,.025) # values of a^2 to test (equivalent to a = 0.32 0.39 0.45 0.50 0.55 0.59)
allpower = seq(.7, .9, .05)

# Make a table to store the N at each power/h^2 combination
powerdf  = data.frame(matrix(NA, nrow=length(allpower), ncol= length(myAA)) )
colnames(powerdf) = myAA
rownames(powerdf) = allpower
method <-2 #ncp method; for search method use method <- 2
j=0
for (i in 1:length(myAA)){

	thisAA = myAA[i]
	print(paste0('AA = ', thisAA))
	for (p in 1:length(allpower)){
	  j<-j+1
		mypower = allpower[p]
		print(paste0('power = ', mypower))
		if(method==1){
		  thispower = power.ACE.test(AA = thisAA, CC = 0, update = c("a_after_dropping_c"), power =mypower, method = "ncp", tryHard = "yes")
		}#
    if(method==2){
      thispower<-power.ACE.test(AA = thisAA, CC = 0, update = c("a_after_dropping_c"), power =mypower, method = "ncp", 
                     tryHard = "yes", search=TRUE)
    }
		thispower$upper<-thisAA
		thispower$lower<-mypower

		powerdf[p, i] = thispower[1] # save the (total) n (pairs)
		if(method==2){
		  w<-which(thispower$power>mypower)[1] #first value greater than specified power
		  lw <- w-1 #previous value
		  #interpolate N to exact power value
		  myfraction <- (mypower-thispower$power[lw])/(thispower$power[w]-thispower$power[lw])
		  interpolatedN <- thispower$N[lw]+myfraction*(thispower$N[w]-thispower$N[lw])
		  powerdf[p,i]<-round(interpolatedN,0)
		}
		# note: an alternative approach is to simulate power empirically
		# thispower = power.ACE.test(AA = thisAA, CC = 0, update = c("a_after_dropping_c"), power =mypower, method = "empirical", tryHard = "yes")
		if (j==1){
		  bigpower <- thispower
		}
		if (j>1){
		  bigpower<-rbind(bigpower,thispower)
		}
	}
}
colnames(bigpower)[3:4]<-c('power','a2')
# The n reported with thispower is total n, MZ and DZ combined
print("N (total pairs of twins) required for a given power (in rows) and heritability (across columns)")
powerdf
write.csv(powerdf, 'powertable.csv')



thispower<-power.ACE.test(AA = thisAA, CC = 0, update = c("a_after_dropping_c"),method = "ncp", 
                          tryHard = "yes", search=TRUE)