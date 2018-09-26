
meanFun <- function(x, group){
  out = tapply(X = x, INDEX = group, FUN = mean)
  return(out)
}

varFun <- function(x, group){
  out = tapply(X = x, INDEX = group, FUN = var)
  return(out)
}

# COMMENT: why it is noRef, when X is given?
ICeDT_fit_noWgt_noRef <- function(Y, X, cellType, fixedCT = NULL, 
                                  fixedCT_rho = NULL, useRho = FALSE, 
                                  borrow4SD = TRUE, maxIter_prop = 100, 
                                  maxIter_PP = 100, RhoConv_CO = 1e-4, 
                                  Subj_CO){
  
  #-----------------------------------------------------#
  # Check Input                                         #
  #-----------------------------------------------------#
  
  if(nrow(Y)!=nrow(X)){
    stop("Expression Inputs do not have the same number of genes!")
  }
  
  if(is.null(fixedCT)){
    message("WARNING: No Fixed Cell Type Present!
            All Proportions will be estimated.")
  } else {
    if(!(fixedCT %in% unique(cellType))){
      stop("Specified fixed cell type label is not
           present in cellType vector!")
    }
    
    if(ncol(Y) != length(fixedCT_rho)){
      stop("Mixture Expressions and Tumor Purity
           mismatch - Different Number of Subjects!")
    } else if(!all(colnames(Y) == names(fixedCT_rho))){
      stop("Mixture Expression and Tumor Purity 
           mismatch! -- Different order of subjects or labels incorrect.")
    }
  }
  
  if(ncol(X) != length(cellType)){
    stop("Number of cell type labels does 
         not match the number of samples in pure references!")
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
  nS  = ncol(Y)
  Qp1 = length(sortCT)
  
  #-----------------------------------------------------#
  # Pure Sample Estimation                              #
  #-----------------------------------------------------#
  # Remnant code which should not influence estimation as
  # seen later. It is not adjusted for the sake of 
  
  # COMMENT: for the sake of what?
  # COMMENT: what if some X is 0?
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
  
  # Edit to ensure tumor contribution is 0
  Z_1[,1] = rep(0, nG)
  
  #----------- Mixture Sample Proportions --------------#
  P_Est = apply(X = rbind(fixedCT_rho,Y), MARGIN = 2, FUN = lmInit,
                Zm=Z_1[,-c(1)], Zt=Z_1[,1])
  Rho_1 = P_Est
  
  #----------- Aberrant Profile Initial ----------------#
  SigmaInit = apply(X = rbind(Y, Rho_1, fixedCT_rho), MARGIN=2,
                    FUN = AbProf_Init, Z = Z_1, nG = nG, Qval = (ncol(Z_1)-1))
  
  Sigma2M_1 = SigmaInit[1,]
  Sigma2A_1 = SigmaInit[2,]
  
  #---- Additional Values ----#
  Pm_1  = rep(0.5,nS)
  
  #-----------------------------------------------------#
  # IMPLEMENT ALGORITHMIC FIT                           #
  #-----------------------------------------------------#
  Rho_0     = Rho_1
  Sigma2M_0 = Sigma2M_1
  Sigma2A_0 = Sigma2A_1
  
  Z_0       = Z_1
  
  Pm_0  = Pm_1

  #------ Update Proportions + -----#
  PropPlus_Out   = PropPlus_Update(Y = Y,Rho_0 = Rho_0,fixedCT_rho = fixedCT_rho,useRho=useRho,
                                   Sigma2M_0 = Sigma2M_0,Sigma2A_0 = Sigma2A_0,
                                   Z_0 = Z_0,Pm_0 = Pm_0,maxIter_PP = maxIter_PP,
                                   maxIter_prop = maxIter_prop,nG = nG,
                                   RhoConv_CO = RhoConv_CO,Subj_CO = Subj_CO)
  
  Rho_1     = PropPlus_Out$Rho
  Sigma2M_1 = PropPlus_Out$Sigma2M
  Sigma2A_1 = PropPlus_Out$Sigma2A
  Pm_1      = PropPlus_Out$Pm
  
  EM_wgt =  HS2_UpdateWgts_All(logY = log(Y),Rho_init = Rho_1,fixed_rho = fixedCT_rho,
                               Sigma2M = Sigma2M_1,Sigma2A = Sigma2A_1,Z = Z_1,p_m = Pm_1)
  
  #-----------------------------------------------------#
  # OUTPUT                                              #
  #-----------------------------------------------------#
  if(useRho){
    fidx = which(sortCT==fixedCT)
    
    IC_Abundance  = cbind(fixedCT_rho,t(Rho_1))
    colnames(IC_Abundance) = c(fixedCT,sortCT[-fidx])
    rownames(IC_Abundance) = colnames(Y)
    
    Fixed_CellType = fixedCT
  } else {
    fidx = which(sortCT==fixedCT)
    
    fixedCT_rho_est = 1-colSums(Rho_1)
    
    IC_Abundance  = cbind(fixedCT_rho_est,t(Rho_1))
    colnames(IC_Abundance) = c(fixedCT,sortCT[-fidx])
    rownames(IC_Abundance) = colnames(Y)
    
    Fixed_CellType = fixedCT
  }
  
  outList = list(IC_Abundance   = IC_Abundance,
                 Fixed_CellType = Fixed_CellType,
                 Sigma2M        = Sigma2M_1,
                 Sigma2A        = Sigma2A_1,
                 Z              = Z_1,
                 CT_VAR         = CT_VAR,
                 P_Consistent   = Pm_1,
                 PP_Consistent  = EM_wgt)
  
  return(outList)
}

PropPlus_Update<- function(Y,Rho_0,fixedCT_rho,useRho,
                           Sigma2M_0,Sigma2A_0,Z_0,Pm_0,
                           maxIter_PP,maxIter_prop,nG,
                           RhoConv_CO,Subj_CO){
  logY = log(Y)
  
  Rho_t1     = Rho_0
  Sigma2M_t1 = Sigma2M_0
  Sigma2A_t1 = Sigma2A_0
  Pm_t1      = Pm_0
  
  for(j in 1:maxIter_PP){
    #----  Reset the Param   ----#
    Rho_t0     = Rho_t1
    Sigma2M_t0 = Sigma2M_t1
    Sigma2A_t0 = Sigma2A_t1
    Pm_t0      = Pm_t1
    
    #---- Update EM Weights  ----#
    EM_wgt = HS2_UpdateWgts_All(logY = logY,Rho_init = Rho_t0,fixed_rho = fixedCT_rho,
                            Sigma2M = Sigma2M_t0,Sigma2A = Sigma2A_t0,Z = Z_0,p_m = Pm_t0)
    
    #---- Update Proportions ----#
    PropP_t1  = HS_UpdatePropn_All(logY = logY,Rho_init = Rho_t0,fixed_rho = fixedCT_rho,
                                   Z = Z_0,maxIter_prop = maxIter_prop,EM_wgt = EM_wgt,
                                   useRho = useRho)
    
    Rho_t1     = PropP_t1$propCurr
    Sigma2M_t1 = PropP_t1$sig2M_Curr
    Sigma2A_t1 = PropP_t1$sig2A_Curr
    
    #---- Update Cons. MarkP ----#
    Pm_t1 = RP_UpdatePM(nG = nG,EM_wgt = EM_wgt)
    
    #---- Check Convergence  ----#
    Rho_Diff = abs(Rho_t0-Rho_t1)
    max_Diff = apply(X = Rho_Diff,MARGIN = 2,FUN = max)
    
    if(sum(max_Diff<RhoConv_CO)>=Subj_CO){
      break
    }
    
    message("Current PropPlus Iter ",j,": Max Diff of ",max(max_Diff))
  }
  
  return(list(Rho = Rho_t1,Sigma2M = Sigma2M_t1,Sigma2A = Sigma2A_t1,Pm = Pm_t1,Iter = j))
}