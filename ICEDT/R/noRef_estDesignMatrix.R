#####################################################################
# Estimating Design Matrix from pure sample expression data         #
#####################################################################
# PROGRAM NAME:                                                     #
#   noRef_estDesignMatrix.R                                         #
# PROGRAMMER:                                                       #
#   Douglas Roy Wilson, Jr.                                         #
# DATE CREATED:                                                     #
#   1/26/2019                                                       #
# LAST EDIT:                                                        #
#   1/26/2019                                                       #
# VERSION:                                                          #
#   R-3.3.1                                                         #
#-------------------------------------------------------------------#
# DESCRIPTION:                                                      #
#   This program is a convenience function designed to separate the #
#   estimation of the design matrix from the ICeD-T modelling       #
#   functions.                                                      #
#####################################################################

meanFun <- function(x, group){
  out = tapply(X = x, INDEX = group, FUN = mean)
  return(out)
}

varFun <- function(x, group){
  out = tapply(X = x, INDEX = group, FUN = var)
  return(out)
}

noRef_estDesign <- function(){
  #RESPONSE: To comment 3, added these checks and edits.
  if(any(is.na(X))){
    stop("Y and X must contain no NA entries. Please remove
         any rows with NA and try again.")
  }
  if(any(Y<0)|any(X<0)){
    stop("As normalized expression values, Y and X should be
         non-negative. Please correct entries less than 0
         and try again.")
  }
  if(any(X<1e-4)){
    message("X or Y contain expression values of 0. Adding
            small correction (1e-5) to ensure log transformation
            is viable.")
    X = X+1e-5
  }
  
  if(ncol(X) != length(cellType)){
    stop("Number of cell type labels does 
         not match the number of samples in pure references!")
  }
  
  if(is.null(fixedCT)){
    fixedCT  = "Dummy_Tumor"
    cellType = c(cellType,"Dummy_Tumor","Dummy_Tumor")
    X        = cbind(X,matrix(0,nrow=nrow(X),ncol=2))
  }
  
  SampTable = table(cellType)
  
  if(min(SampTable)<2){
    stop("At least two pure samples are needed 
         for each cell type.")
  }
  
  #-----------------------------------------------------#
  # EXTRACTING SAMPLE INFO                              #
  #-----------------------------------------------------#
  # Assumes some pure tumor subjects! If this is not met,
  # then the program is stopped above and user is counseled.
  
  nnct = tapply(cellType,cellType,length)
  ntot = sum(nnct)
  
  sortCT = sort(unique(cellType))
  
  nG  = nrow(X)
  nP  = ncol(X)
  Qp1 = length(sortCT)
  
  logX = log(X) 
  
  CT_MU  = t(apply(X = logX, MARGIN = 1, FUN = meanFun, 
                   group = cellType))
  
  CT_VAR = t(apply(X = logX, MARGIN = 1, FUN = varFun,
                   group = cellType))
  
  
  #------------------ Quick Check ----------------------#
  
  if(any(colnames(CT_MU)!=sortCT)){
    message("Problems with tapply order!")
  }
  
  if(any(colnames(CT_VAR)!=sortCT)){
    message("Problems with tapply order!")
  }
  
  if(any(names(nnct)!=sortCT)){
    message("Problems with tapply order!")
  }
  
  #----------------- Borrow For SD ---------------------#
  if(borrow4SD){
    X1 = logX
    ntot_p = length(cellType) - length(which(cellType==fixedCT))
    
    for(i in 1:Qp1){
      if(sortCT[i] == fixedCT){ next }
      wwi = which(cellType == sortCT[i])
      X1[,wwi] = logX[,wwi] - CT_MU[,i]
    }
    
    varXPop = apply(X1[,-which(cellType==fixedCT)], 1, var)
    
    for(i in 1:Qp1){
      if(sortCT[i] == fixedCT){ next }
      wi = (nnct[i])/(ntot_p)
      CT_VAR[,i] = wi*CT_VAR[,i] + (1-wi)*varXPop
    }
  }
  
  #----------------- Correct Order ---------------------#
  # Rearrange order so that fixed cell type (tumor) is the first one
  fIdx   = which(colnames(CT_MU)==fixedCT)
  CT_MU  = cbind(CT_MU[,fIdx],CT_MU[,-c(fIdx)])
  CT_VAR = cbind(CT_VAR[,fIdx],CT_VAR[,-c(fIdx)])
  colnames(CT_MU)[1]  = fixedCT
  colnames(CT_VAR)[1] = fixedCT
  
  #-----------------------------------------------------#
  # INITIALIZATION                                      #
  #-----------------------------------------------------#
  CT_VAR_1 = CT_VAR
  Z_1      = exp(CT_MU+CT_VAR/2)
  
  # Edit to ensure tumor (fixed cell type) contribution is 0
  Z_1[,1] = rep(0, nG)
  
  # Return the estimated reference matrix with the column 
  # pertaining to the tumor cell type set to 0. 
  return(Z_1[,-c(1)])
}