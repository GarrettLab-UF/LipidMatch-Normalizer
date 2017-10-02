#########################################################
# LipidMatch Quant                                      #
# Author: Jason Cochran                                 #
# Chemistry/moral support: Jeremy Koelmel               #
#########################################################

rm( list = ls() )

#### GUI starts here ###############################################################################
if("tcltk2" %in% rownames(installed.packages()) == FALSE) {install.packages("tcltk2")}
library(tcltk2)

wd <- tk_choose.dir(default = "", caption = "Select directory")
setwd(wd)
done <- tclVar(0)

numAdducts <- NULL
numValues <- NULL
rt_tolerance <- NULL
mz_tolerance <- NULL
featureTable_loc <- NULL
intStd_loc <- NULL
output <- NULL

RTCol <- NULL
mzCol <- NULL
sampleStartCol <- NULL
sampleEndCol <- NULL
classIDCol <- NULL
adductIDCol <- NULL
numericDataStart_row <- NULL
sampleGrouping_row <- NULL
normalizeByWeights <- NULL

weights_col <-NULL

GUILauncher <- function() {
  
  output <<- tclVar("Default Output Folder")
  mz_tolerance <<- tclVar(".005")
  rt_tolerance <<- tclVar(".15")
  normalizeByWeights <<- tclVar("0")
  
  numAdducts <<- tclVar("5")
  numValues <<- tclVar("18")
  sampleStartCol <<- tclVar("7")
  sampleEndCol <<- tclVar("20")
  mzCol <<- tclVar("3")
  RTCol <<- tclVar("2")
  classIDCol <<- tclVar("21")
  adductIDCol <<- tclVar("22")
  numericDataStart_row <<- tclVar("3")
  sampleGrouping_row <<- tclVar("2")

  tt <- tktoplevel()
  tkwm.title(tt,"LipidMatchQuant Settings")
  numValues.entry <- tkentry(tt, textvariable= numValues)
  output.entry <- tkentry(tt, textvariable = output)
  mz_tolerance.entry <- tkentry(tt, textvariable= mz_tolerance)
  rt_tolerance.entry <- tkentry(tt, textvariable= rt_tolerance)
  normalizeByWeights.entry <- tk2checkbutton(tt, text = "Normalize data (e.g. by weight or protein)")
  
  numAdducts.entry <- tkentry(tt, textvariable= numAdducts)
  sampleStartCol.entry <- tkentry(tt, textvariable  = sampleStartCol)
  sampleEndCol.entry <- tkentry(tt, textvariable = sampleEndCol)
  mzCol.entry <- tkentry(tt, textvariable = mzCol)
  RTCol.entry <- tkentry(tt, textvariable = RTCol)
  classIDCol.entry <- tkentry(tt, textvariable = classIDCol)
  adductIDCol.entry <- tkentry(tt, textvariable = adductIDCol)
  numericDataStart_row.entry <- tkentry(tt, textvariable = numericDataStart_row)
  sampleGrouping_row.entry <- tkentry(tt, textvariable = sampleGrouping_row)

  submit <- function() {
    output <<- tclvalue(output)
    mz_tolerance <<- as.numeric(tclvalue(mz_tolerance))
    rt_tolerance <<- as.numeric(tclvalue(rt_tolerance))
    normalizeByWeights <<- as.character(tclvalue(normalizeByWeights))
    numAdducts <<- as.numeric(tclvalue(numAdducts))
    numValues <<- as.numeric(tclvalue(numValues))
    sampleStartCol <<- as.numeric(tclvalue(sampleStartCol))
    sampleEndCol <<- as.numeric(tclvalue(sampleEndCol))
    mzCol <<- as.numeric(tclvalue(mzCol))
    RTCol <<- as.numeric(tclvalue(RTCol))
    classIDCol <<- as.numeric(tclvalue(classIDCol))
    adductIDCol <<- as.numeric(tclvalue(adductIDCol))
    numericDataStart_row <<- as.numeric(tclvalue(numericDataStart_row))
    sampleGrouping_row <<- as.numeric(tclvalue(sampleGrouping_row))
    tclvalue(done) <- 1
    tkdestroy(tt)
  }
  submit.but <- tkbutton(tt, text="Run Quantification", command=submit)
  
  FT_loc <- function() {
    featureTable_loc <<-tk_choose.files(default = "", caption = "Select Feature Table", multi = FALSE, filters = NULL, index = 1)
  }
  FT_loc.but <- tkbutton(tt, text="Select Feature Table", command = FT_loc)
  
  IS_loc <- function() {
    intStd_loc <<- tk_choose.files(default = "", caption = "Select Internal Standard Table", multi = FALSE, filters = NULL, index = 1)
  }
  IS_loc.but <- tkbutton(tt, text="Select Internal Standard Table", command = IS_loc)
  
  tkgrid(tklabel(tt, text="Select Feature Table with button") , FT_loc.but)
  tkgrid(tklabel(tt, text="Select IS Table with button") ,IS_loc.but)
  tkgrid(tklabel(tt, text=""))
  tkgrid(tklabel(tt, text="Output folder name (created automatically): "), output.entry)
  tkgrid(tklabel(tt,text="m/z Tolerance:"), mz_tolerance.entry)
  tkgrid(tklabel(tt,text="RT Tolerance:"), rt_tolerance.entry)
  tkgrid(tklabel(tt, text=""))
  tkgrid(tklabel(tt,text="Number of Adducts:"), numAdducts.entry)
  tkgrid(tklabel(tt, text="Sample Start Column:"), sampleStartCol.entry)
  tkgrid(tklabel(tt, text="Sample End Column:"), sampleEndCol.entry)
  tkgrid(tklabel(tt, text="m/z Column:"), mzCol.entry)
  tkgrid(tklabel(tt, text="RT Column:"), RTCol.entry)
  tkgrid(tklabel(tt, text="Class ID Column:"), classIDCol.entry)
  tkgrid(tklabel(tt, text="Adduct ID Column:"), adductIDCol.entry)
  tkgrid(tklabel(tt, text="Numeric Data Start Row:"), numericDataStart_row.entry)
  tkgrid(tklabel(tt, text="Sample Grouping Row:"), sampleGrouping_row.entry)
  tkgrid(tklabel(tt, text=""))
  tkgrid(submit.but, pady = 8, columnspan=2)
}

GUILauncher()
tkwait.variable(done)
if( !dir.exists(output)) {
  dir.create(output)
}

##### GUI ends here ###########################################################################

# Data import and cleaning
InternalStandard <- as.matrix( read.csv(basename(intStd_loc), header = F, stringsAsFactors = F) )
# Weights will be used later...
weights <- InternalStandard[2,]
InternalStandard <- InternalStandard[-2,]
titles <- as.matrix(read.csv(basename(intStd_loc), header = F, stringsAsFactors = F))[2,3:(2+numAdducts)]
InternalStandard[1,3:(2+numAdducts)] <- titles
titles <- InternalStandard[1,]
InternalStandard <- as.data.frame(InternalStandard[-1,], col.names = names(titles), stringsAsFactors = F )
colnames(InternalStandard) <- titles

FeatureTable <- read.csv(basename(featureTable_loc), header = T, stringsAsFactors = F)
grouping <- FeatureTable[sampleGrouping_row-1,]

for(i in 2:(numericDataStart_row - 1) ) {
  FeatureTable <- FeatureTable[-1,]
}


row <- 1

toleranceCheck = function(target, givenValue, tolerance) {
  if( ((target + tolerance) > givenValue) && ((target - tolerance) < givenValue) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Redoing locating the standards in the feature table
matches <- matrix(data=NA, nrow=0, ncol = (ncol(FeatureTable) + 2))

FindStandards = function(selISRow) {
  # For each row of Internal Standard table search for anything within the RT tolerance window
  potentialStandards <- subset(FeatureTable, FeatureTable[,RTCol] <= ( as.numeric(selISRow[2]) + rt_tolerance) )
  potentialStandards <- subset(potentialStandards, potentialStandards[,RTCol] >= (as.numeric(selISRow[2]) - rt_tolerance) )
  # Once we have that list... we will then search it multiple times (depending on the # of adducts) using the m/z of each adduct
  for(i in 1:numAdducts) {
    # If a match is found for an adduct, we will label it in the matches table. If there are multiple, choose the one with the closest RT.
    potentialAdductMatches <- subset(potentialStandards, potentialStandards[,mzCol] <= (as.numeric(as.numeric(selISRow[2+i])) + mz_tolerance) )
    potentialAdductMatches <- subset(potentialAdductMatches, potentialAdductMatches[,mzCol] >= (as.numeric(selISRow[2+i]) - mz_tolerance))
    matchAdduct <- names(selISRow[2+i])
    if( nrow(potentialAdductMatches) >= 2) { # We only need 1 IS... remove all others
      # Reduce the possibilities using the closest RT of the matches to the original IS
      closestRow <- which( abs( as.numeric(potentialAdductMatches[,RTCol]) - as.numeric(selISRow[2]) ) == min(abs( as.numeric(potentialAdductMatches[,RTCol]) - as.numeric(selISRow[2]) )) )
      potentialAdductMatches <- potentialAdductMatches[closestRow,]
    }
    # Now we should check if we actually have any sort of match
    if( nrow(potentialAdductMatches) == 1 ) {
      # Need to save the rows of the feature table that are matched... add columns on to them for "matchClass" and "matchAdduct"
      FeatureTable <<- FeatureTable[ which(FeatureTable[,1] != potentialAdductMatches[,1]) ,]
      matchClass <- selISRow[1]
      currentMatch <- cbind(potentialAdductMatches, matchClass, matchAdduct)
      matches <<- rbind(matches, currentMatch)
    } else if( nrow(potentialAdductMatches) == 0) {
      # Do nothing.
    } else {
      print("Multiple features identified as matches for a single internal standard. Please contact the developers.")
    }
  }
}

invisible( apply(InternalStandard, 1, FindStandards) )

rm(mz_tolerance, row, rt_tolerance, titles, FindStandards)

#############################################
############### PART 2 & 3 ##################
#############################################

# Make a table of just our matched Intd Standards to make quantifying easier
write.table(matches, file = paste(output ,"standardsFound.csv", sep = "/"), sep = ",", col.names = TRUE, row.names = FALSE)

# Setup a dataframe to store all the classes we need
quantClasses <- as.list( strsplit( as.character( InternalStandard$Classes), split = " " ) )
quantClasses <- sapply( quantClasses , '[', seq(max(sapply(quantClasses, length))))
if(!is.matrix(quantClasses)){
  quantClasses <- as.matrix(quantClasses)
} else {
  quantClasses <- t(quantClasses)
}
numQuantClasses <- ncol(quantClasses)
classes <- cbind.data.frame(InternalStandard[,1:(2+numAdducts)], InternalStandard[,(numAdducts+3):length(InternalStandard)], quantClasses )

isQuantified <- matrix(data=0, nrow = nrow(FeatureTable), ncol = 1)
IS_Species <- matrix(data = NA, nrow = nrow(FeatureTable), ncol = 1)
IS_Adduct <- matrix(data = NA, nrow = nrow(FeatureTable), ncol = 1)
quantifiedAmounts <- cbind(FeatureTable, isQuantified, IS_Species, IS_Adduct)
rm(quantClasses, IS_Species, IS_Adduct)

quantifier_actual = function(rawInput, curStandard, curIntStd) {
  rowNum <- which(quantifiedAmounts[,1] == as.numeric(rawInput[1]) )
  quantifiedAmounts$isQuantified[rowNum] <<- 1
  quantifiedAmounts$IS_Species[rowNum] <<- curIntStd[[1]]
  quantifiedAmounts$IS_Adduct[rowNum] <<- as.character(unlist(curStandard$matchAdduct))
  i <- sampleStartCol
  j <- 1
  while(i <= sampleEndCol) {
    quantifiedAmounts[rowNum, i] <<- ( as.double(curIntStd[3+numAdducts+j]) / as.double(curStandard[i]) ) *  as.double(rawInput[i])
    # For ^: (IS amount / IS Signal ) * Sample signal
    i <- i + 1
    j <- j + 1
  }
}

quantifier_multi_IS = function(class_sel) {
  masterClassList <- NULL
  # Select anything with the class we want
  masterClassList <- subset(FeatureTable, FeatureTable[,classIDCol] == class_sel[[1]])
  # Narrow down by the adduct, use a narrowed down adduct IS list for this purpose
  #     For each adduct we have of that class...
  #         Find all standards which are available for that adduct (check matches table for this data)
  #         For all the standards available make a table of what RT goes with each standard and its respective adduct
  for(i in 1:numAdducts) {
    curAdduct <- colnames(InternalStandard)[[2+i]]
    curStandards <- matrix(data = NA, nrow = 0, ncol = 2)
    currentFeatureList <- subset(masterClassList, masterClassList[,adductIDCol] == curAdduct)
    noStandardFound <- TRUE
    if( nrow(currentFeatureList) != 0 ) {
      # class_sel[2] gives a list of all standards that we can actually use for this CLASS... now we need to find out what we can use for this adduct
      for(numIS in 1:length(class_sel[[2]])) {
        # For each IS candidate we need to find out if it has the adduct we want is available
        tempStandards <- subset(matches, matches$matchClass == class_sel[[2]][[numIS]])
        tempStandards <- subset(tempStandards, tempStandards$matchAdduct == curAdduct )
        if( nrow(tempStandards) == 1 ) {
          temp_insert_row <- matrix(data = NA, nrow = 1, ncol = 2)
          temp_insert_row[,1] <- as.character(unlist(tempStandards$matchClass[[1]]))
          temp_insert_row[,2] <- tempStandards[,4]
          curStandards <- rbind(curStandards, temp_insert_row)
          noStandardFound <- FALSE
        } else {
          # print("Error: More than one standard has been detected when checking for adduct availability.")
        }
      }
      # Have a list of all standards that quantify this particular adduct and class
      # Determine what retention time ranges are associated with what standard available for that class (based on the adduct)
      # Order the curStandards table by RT...
      if(nrow(curStandards) != 1 ) {
        curStandards <- curStandards[,order( as.numeric(curStandards[,2]) )]
      }
      # browser()
      RT_cutoffs <- matrix(data = NA, nrow = 1, ncol = ( nrow(curStandards) + 1) )
      RT_cutoffs[ncol(RT_cutoffs)] = 10000000000000
      RT_cutoffs[1] = 0
      #First we create the RT cutoff table
      if( ncol(RT_cutoffs) > 2) {
        for( RT_value in 1:( nrow(curStandards) - 1 ) ) {
          # RT_cutoffs[1+RT_value] <<- ( as.numeric( curStandards[RT_value,2]  ) + as.numeric( curStandards[RT_value+1,2] ) ) / 2
          RT_cutoffs[1+RT_value] <- ( as.numeric( curStandards[RT_value,2]  ) + as.numeric( curStandards[RT_value+1,2] ) ) / 2
        }
      }
      # Then we quantify
      #   First we take everything that is the correct adduct and class
      #   Then we select things that fall within each range and send those things with that standard to be quantified
      if( !noStandardFound ) {
        for(standard in 1:(nrow(curStandards)) ) {
          # First select our standard from matches
          curStandard <- subset(matches, matches$matchClass == curStandards[standard,1] )
          curStandard <- subset(curStandard, curStandard$matchAdduct == curAdduct)
          # Now select the internal standard from the IS table
          curIntStd <- subset(InternalStandard, InternalStandard[,1] == curStandards[standard,1])
          # Then we will get all of the features we can quant with this standard
          selected_subset <- subset(currentFeatureList, currentFeatureList[,RTCol] <= RT_cutoffs[standard+1] )
          selected_subset <- subset(selected_subset, selected_subset[,RTCol] >= RT_cutoffs[standard] )
          apply( selected_subset, 1, quantifier_actual, curStandard = curStandard, curIntStd = curIntStd)
        }
      }
    }
  }
}

# Make a list of classes that you can use to quantify each
#     First identify all unique lipid classes names... then identiy using row numbers what standards can be used for each
# column start reference 2+numAdducts+( sampleEndCol - sampleStartCol + 1 ) + i
# is iteration number
class_standards <- data.frame(Class = as.character(), IS_list = I(as.list(as.character())), stringsAsFactors = FALSE )

for(i in (2 + numAdducts + (sampleEndCol - sampleStartCol + 1) + 2):ncol(classes) )  {
  for( j in 1:(nrow(classes)) ) {
    if( !is.na(classes[j,i])) {
      count <- sum(class_standards[,1] == classes[j,i])
      if( count == 0 ) {
        temp_new <- data.frame(Class = as.character(), IS_list = I(as.list(as.character())), stringsAsFactors = FALSE )
        temp_new[1,] <- classes[j,i] # class abbreviation eg. CER, LS
        temp_new[1,2] <- list(classes[j,1]) # IS name
        class_standards <<- rbind(class_standards, temp_new)
      } else { # Standard already accounted for. Add another possible IS for use
        rowNum <- which(class_standards[,1] == classes[j,i])
        class_standards[[rowNum, 2]] <- c( class_standards[rowNum,2][[1]], classes[j,1] )
      } 
    }
  }
}

apply( class_standards, 1, quantifier_multi_IS)

temp_init <- subset(quantifiedAmounts, quantifiedAmounts$isQuantified == 1)
rm( quantifier_actual)

#############################################
# Quantify anything that we didn't find a int std for but that we identified
############### PART 4 ##################
  
unQuantified <- quantifiedAmounts[ quantifiedAmounts$isQuantified == 0, ]
unQuantified[,classIDCol][unQuantified[,classIDCol] == ""] <- NA
unQuantified[,adductIDCol][unQuantified[,adductIDCol] == ""] <- NA
unQuantified <- unQuantified[complete.cases(unQuantified[,classIDCol]),]
unQuantified <- unQuantified[complete.cases(unQuantified[,adductIDCol]),]

# tempQuantifiers <- matrix( unique( unQuantified$Class_At_Max_Intensity ) )
tempQuantifiers <- matrix( unique( unQuantified[,classIDCol] ) )
desiredQuantifiers <- matrix(data = NA, nrow = 0, ncol = 2)

# Now lets get the unique Adducts for each class as well
i <- 1
while(i <= nrow(tempQuantifiers) ) {
  tempSubset <- subset(unQuantified, unQuantified[,classIDCol] == tempQuantifiers[i] )
  tempAdducts <- matrix( unique( tempSubset[,adductIDCol] ))
  j <- 1
  while(j <= nrow(tempAdducts)) {
    tempRow <- matrix(data=NA, nrow = 1, ncol = 2 )
    tempRow[1,1] <- tempQuantifiers[i]
    tempRow[1,2] <- tempAdducts[j]
    desiredQuantifiers <- rbind(desiredQuantifiers, tempRow )
    j <- j + 1
  }
  i = i + 1
}

desiredQuantifiers

rm(tempQuantifiers, j, i, tempSubset, tempAdducts, tempRow)

# numValues - number of values in columns to quantify
# numAdduct - number of adducts in IntStd table
# classes - subseted and cleaned view of internal standards

quantifier = function(rawInput, sel_IS, curStandard, score) {
  rowNum <- which(quantifiedAmounts[,1] == as.numeric(rawInput[1]) )
  quantifiedAmounts$isQuantified[rowNum] <<- score
  quantifiedAmounts$IS_Species[rowNum] <<- sel_IS[[1]]
  quantifiedAmounts$IS_Adduct[rowNum] <<- as.character(unlist(curStandard$matchAdduct))
  
  i <- sampleStartCol
  j <- 1
  while(i <= sampleEndCol) {
    # quantifiedAmounts[rowNum, i] <<- ( as.double(curIntStd[3+numAdducts+j]) / as.double(curStandard[i]) ) *  as.double(rawInput[i])
    # Below is old line from this section; above is line from score 1 quantification
    # quantifiedAmounts[rowNum, i] <<- ( as.double(sel_IS[3+numAdducts+j-1]) / as.double(curStandard[i]) ) *  as.double(rawInput[i])
    quantifiedAmounts[rowNum, i] <<- ( as.double(sel_IS[3+numAdducts+j]) / as.double(curStandard[i]) ) *  as.double(rawInput[i])
    # For ^ (IS amount / IS Signal ) * Sample signal
    i <- i + 1
    j <- j + 1
  }
}

comparator = function(sel_group) {
  subset_sm <- subset(unQuantified, unQuantified[,classIDCol] == sel_group[1] )
  subset_sm <- subset(subset_sm, subset_sm[,adductIDCol] == sel_group[2] )
  if( nrow(subset_sm) >= 1) {
    avgRT <- mean( subset_sm[,RTCol] )
    # Sort out all classes that don't have value for the adduct we want... if we get list of length(0) then
    # just get the closest RT match
    adduct_desired <- which(names(classes[,3:(3+numAdducts-1)]) == sel_group[2] )
    subset_classes <- classes[ !is.na(classes[,adduct_desired + 2]), ]
    if( nrow(classes[ !is.na(classes[,adduct_desired + 2]), ]) >= 1 ){
      # Now find standard with closest value
      possibleISSubset <- subset(matches, matches$matchAdduct == sel_group[[2]])
      avgIntStd <- which( abs( as.numeric(possibleISSubset[,RTCol])-avgRT) == min(abs( as.numeric(possibleISSubset[,RTCol])-avgRT)) )
      # if nothing gets selected in the above stage then error out
      if( length(avgIntStd) != 0) {
        # Set avgIntStd to be the row from IntStd.csv we will quantify with
        avgIntStd <- possibleISSubset[avgIntStd,]
      } else {
        print("Error: Found 0 Int Std for this class")
      }
      # Now let's use the standard we know we have to find the standard inside the IS table...
      curStandard <- avgIntStd
      avgIntStd <- which( subset_classes[,1] == curStandard$matchClass)
      avgIntStd <- subset_classes[avgIntStd,]
      
      # Now let's quantify according to that standard; quantifiedAmounts will still be our target...
      if( nrow(curStandard) != 0 ) {
        apply(subset_sm, 1, quantifier, sel_IS = avgIntStd, curStandard = curStandard, score = 2 )
      }
    } else {
      subset_classes <- classes
      avgIntStd <- which( abs( as.numeric(subset_classes$RT)-avgRT) == min(abs( as.numeric(subset_classes$RT)-avgRT)) )
      avgIntStd <- subset_classes[avgIntStd,]
      curStandard <- which( matches$matchClass == avgIntStd[,1]  )
      curStandard <- matches[curStandard,]
      # Need to refine selection of curStandard to account for m/z average
      mzAvg <- mean(subset_sm[,mzCol])
      mzAvgFeature <- which( abs( as.numeric(curStandard[,mzCol])-mzAvg) == min(abs( as.numeric(curStandard[,mzCol])-mzAvg)) )
      curStandard <- curStandard[mzAvgFeature,]
      apply(subset_sm, 1, quantifier, sel_IS = avgIntStd, curStandard = curStandard, score = 3 )
    }
  } else {
    # print("No matches identified lipids were found to quantify using our current standard.")
  }
}

apply(desiredQuantifiers, 1, comparator)

temprow <- matrix(c(rep.int(NA,length(quantifiedAmounts))), nrow=1, ncol=length(quantifiedAmounts))
experimentalDesign <- data.frame(temprow)
colnames(experimentalDesign) <- colnames(quantifiedAmounts)
for(i in sampleStartCol:sampleEndCol) {
  experimentalDesign[,i] = grouping[,i]
}

quantifiedAmounts <- rbind(experimentalDesign, quantifiedAmounts)

write.table(quantifiedAmounts, file = paste(output, paste(substr(basename(featureTable_loc),1,nchar(basename(featureTable_loc))-4),"Quant.csv", sep = "_"), sep = "/"), sep = ",", col.names = TRUE, row.names = FALSE)
temp <- subset(quantifiedAmounts, quantifiedAmounts$isQuantified == 2 || quantifiedAmounts$isQuantified == 3)
print(paste("Quantified with score of 1: ", nrow(temp_init), sep = "" ) )
print(paste("Quantified with score of 2 or 3: ", nrow(temp), sep = "" ) )

grouping <- quantifiedAmounts[1,]
quantifiedAmounts <- quantifiedAmounts[-1,]

print("Quantification is complete")
