# Shelli Kesler 6/5/19
# Connectome-based predictive modeling with split-sample CV
# Adapted from Matlab code Copyright 2015 Xilin Shen and Emily Finn
# Reference: Finn ES, et al. Nature Neuroscience 18, 1664-1671.

# Uses edge sums for each of 8 networks (Shen, X, et al. 2013. NeuroImage, 82, pp. 403-15.)
# Conducts simple linear regression as well as random forest regression, support vector regression
# with a linear kernel and support vector regression with an RBF kernel 

# NOTE: fdata (Matlab nxnxm connectivity matrix) var name must be Z
# bdata = behavioral data: assumes this is preloaded in the R environment
# p = proportion of sample to use for training set
# Usage e.g. CPM('Z57.mat',cog,.75)

CPM50 <- function(fdata,bdata,p){
  library(R.matlab)
  library(randomForest)
  library(e1071)
  
  X = readMat(fdata)
  bdata = as.matrix(bdata)
  nsubs=dim(X$Z)[3]
  nnodes=dim(X$Z)[2]
  
  thresh = .05
  
  nx = 1:nsubs
  inTrain = sample(nx,round(p*length(nx)))
  inTest = nx[-inTrain]
  
  all_mats = X$Z
  
  train_mats = array(NA, dim=c(nnodes,nnodes,length(inTrain)))
  
  for (t in 1:length(inTrain)){
    train_mats[,,t] = all_mats[,,inTrain[t]]
  }
  
  train_behav = bdata[inTrain]
  
  test_mats = array(NA, dim=c(nnodes,nnodes,nsubs-length(inTrain)))
  
  for (ts in 1:length(inTest)){
    test_mats[,,ts] = all_mats[,,inTest[ts]]
  }
  
  test_behav = bdata[inTest]
  
  predictedValues = matrix(data = NA, nrow=length(inTest),ncol=4)
  netLabs <- c('Net1','Net2','Net3','Net4','Net5','Net6','Net7','Net8')
  
  train_vcts = matrix(data=NA,ncol=length(inTrain),nrow=nnodes^2)
   
      for (s in 1:dim(train_mats)[3]){
        train_vcts[,s] = as.vector(train_mats[,,s])
      }
    
    # correlate edges with behavior
    r_mat = cor(t(train_vcts),train_behav) #rcorr, corr.test, etc. don't work with NAs
    r2 = r_mat^2
    p_mat = 1-pbeta(r2,1/2,(length(inTrain)-3)/2) #gives same pvals as Matlab
    
    r_mat = matrix(r_mat,nrow=nnodes,ncol=nnodes)
    p_mat = matrix(p_mat,nrow=nnodes,ncol=nnodes)

    # set threshold and define masks
    pos_mask = matrix(data=0,nrow = nnodes, ncol = nnodes) 
    pos_mask[(p_mat < thresh)] = 1
    
    # get sum of all edges for each network in TRAIN subs (divide by 2 to control for the
    # fact that matrices are symmetric)
    
    posEdge = array(NA, dim=c(nnodes,nnodes,length(inTrain)))
    net_sumpos = matrix(nrow=length(inTrain),ncol=8)
    
      for (ss in 1:(length(inTrain))){   
        posEdge[,,ss] = train_mats[,,ss]*pos_mask
        
        net_sumpos[ss,1] = sum(posEdge[1:29,1:29,ss],na.rm = TRUE)/2
        net_sumpos[ss,2] = sum(posEdge[30:64,30:64,ss],na.rm = TRUE)/2
        net_sumpos[ss,3] = sum(posEdge[65:83,65:83,ss],na.rm = TRUE)/2
        net_sumpos[ss,4] = sum(posEdge[84:173,84:173,ss],na.rm = TRUE)/2
        net_sumpos[ss,5] = sum(posEdge[174:223,174:223,ss],na.rm = TRUE)/2
        net_sumpos[ss,6] = sum(posEdge[224:241,224:241,ss],na.rm = TRUE)/2
        net_sumpos[ss,7] = sum(posEdge[242:250,242:250,ss],na.rm = TRUE)/2
        net_sumpos[ss,8] = sum(posEdge[251:268,251:268,ss],na.rm = TRUE)/2
      }
    
    colnames(net_sumpos) <- netLabs
    # remove any zero sum columns
    inx1 = which(apply(net_sumpos, 2, sum) == 0)
    
    if(length(inx1) != 0){
      net_sumpos = net_sumpos[, -inx1]
    } 
    
    
    # TEST data
    testnet_sumpos = matrix(nrow=length(inTest),ncol=8)
    TestposEdge = array(NA, dim=c(nnodes,nnodes,length(inTest)))
    
    for (sss in 1:(length(inTest))){   
      TestposEdge[,,sss] = test_mats[,,sss]*pos_mask
      
      testnet_sumpos[sss,1] = sum(TestposEdge[1:29,1:29,sss],na.rm = TRUE)/2
      testnet_sumpos[sss,2] = sum(TestposEdge[30:64,30:64,sss],na.rm = TRUE)/2
      testnet_sumpos[sss,3] = sum(TestposEdge[65:83,65:83,sss],na.rm = TRUE)/2
      testnet_sumpos[sss,4] = sum(TestposEdge[84:173,84:173,sss],na.rm = TRUE)/2
      testnet_sumpos[sss,5] = sum(TestposEdge[174:223,174:223,sss],na.rm = TRUE)/2
      testnet_sumpos[sss,6] = sum(TestposEdge[224:241,224:241,sss],na.rm = TRUE)/2
      testnet_sumpos[sss,7] = sum(TestposEdge[242:250,242:250,sss],na.rm = TRUE)/2
      testnet_sumpos[sss,8] = sum(TestposEdge[251:268,251:268,sss],na.rm = TRUE)/2
    }
    
    colnames(testnet_sumpos) <- netLabs
    
    if(length(inx1) != 0){
      testnet_sumpos = testnet_sumpos[, -inx1]
    } 
    
    # Predict scores in TRAIN and apply to TEST
    rfFit <- randomForest(net_sumpos,train_behav,mtry=3,ntree=500,replace=TRUE)
    predictedValues[,1] <- predict(rfFit,testnet_sumpos)
    tempTrain <- as.data.frame(cbind(net_sumpos,train_behav))
    lrFit <- lm(train_behav ~ .,data=tempTrain)
    tempTest <- as.data.frame(cbind(testnet_sumpos,test_behav))
    names(tempTest)[9] <- "train_behav"
    predictedValues[,2] <- predict(lrFit,tempTest)
    svr.lnFit <- svm(train_behav ~ ., tempTrain, kernel = "linear")
    predictedValues[,3] <- predict(svr.lnFit,tempTest)
    svr.rbfFit <- svm(train_behav ~ ., tempTrain, kernel = "radial")
    predictedValues[,4] <- predict(svr.rbfFit,tempTest)
  
  # Evaluate model performance
  RP <- matrix(data=NA,ncol=2,nrow=4)
  R1 <- cor.test(predictedValues[,1],test_behav)
  R2 <- cor.test(predictedValues[,2],test_behav)
  R3 <- cor.test(predictedValues[,3],test_behav)
  R4 <- cor.test(predictedValues[,4],test_behav)
  RP[1,1] <- R1$estimate; RP[1,2] <- R1$p.value
  RP[2,1] <- R2$estimate; RP[2,2] <- R2$p.value
  RP[3,1] <- R3$estimate; RP[3,2] <- R3$p.value
  RP[4,1] <- R4$estimate; RP[4,2] <- R4$p.value
  
  modelLabs <- c("RF","LR","SVRln","SVRrbf")
  colnames(bdata) <- "bdata"
  colnames(predictedValues) <- modelLabs
  predictedVals <<- cbind(predictedValues,test_behav)
  
  rownames(RP) <- modelLabs
  colnames(RP) <- c("R","pval")
  RP <<- RP
  return(RP)
}#end function
