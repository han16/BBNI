#################################
library(bitops)
library(stats)
start.time=date()
ptm <- proc.time()
##################################################################################
update.ancestor_matrix<-function(incid_matrix) #  for a given incidence matrix, generate corresponding ancestor matrix
{
 ances_matrix<-matrix(0, nrow=nrow(incid_matrix), ncol=ncol(incid_matrix))  # important! every time before updating, ancestor matrix should be cleared.
 ances_matrix<-incid_matrix
 for ( ii in 1: nrow(ances_matrix))
 {
  for (i in 1: nrow(ances_matrix))
    for (j in 1:ncol(ances_matrix))
     for (k in 1: nrow(ances_matrix))
       if (ances_matrix[i,j]==1 & ances_matrix[j,k]==1)      # it should be updated num.node times
          ances_matrix[i,k]<-1
 }
 return(ances_matrix)
}
#########################################################################################
check.ances.matrix<-function(ances_matrix) # check whether there are loops in the whole network by checking ancestor matrix
 {
   loop<-0
   for (i in 1:nrow(ances_matrix))
    if (ances_matrix[i,i]==1)
      loop<-loop+1
   return(loop)
 }
#########################################################################################
#########################################################################################
Prop_Trans_Func_Matrix<-function(prop_incid_matrix)      # based on incidence matrix, define transition function matrix
  {
       prop_trans_func_matrix<-prop_incid_matrix; jj<-numeric()
             for (i in 1:nrow(prop_incid_matrix))
               {
                 if (sum(prop_incid_matrix[i,])==1)    #pairwise relation
                   for (j in 1: nrow(prop_incid_matrix))
                     if (prop_incid_matrix[i,j]==1)
                       prop_trans_func_matrix[i,j]<-10+sample.int(2,1,replace=FALSE)

                 if (sum(prop_incid_matrix[i,])==2)    # triplet relation
                  {     j1<-1
                   for (j in 1:ncol(prop_incid_matrix))
                     if (prop_incid_matrix[i,j]==1)
                       {
                        jj[j1]<-j
                        j1<-j1+1
                       }
                    prop_trans_func_matrix[i,jj[1]]<-sample.int(10,1,replace=FALSE)
                    prop_trans_func_matrix[i,jj[2]]<-prop_trans_func_matrix[i,jj[1]]
                   }
                 }
      return(prop_trans_func_matrix)
  }
#################################################################################
#################################################################################
Error_LLH<-function(TRFUM)   # compute error-likelihood, this function depends on the TRansition FUnction Matrix
   {                                             # error_prior is a vector not a matrix
        InOutPair<-list(); ii<-0; root_node<-numeric(); root<-0
        for ( i in 1:nrow(TRFUM))          #  i row   find each input output combination
         {
          if (sum(TRFUM[i,])>0)
            {
             ii<-ii+1
             p<-numeric(); k<-1; jj<-numeric()
             for (j in 1:nrow(TRFUM))  # j column
               if (TRFUM[i,j]!=0)
                 {
                  p[k]<-TRFUM[i,j];jj[k]<-j
                  InOutPair[[ii]]<-c(i,jj,p[k])      # note that InOutPair is a list with vector as its entries.
                  k<-k+1                             # this determines the directionality, i.e. g_k=f(g_i,g_j), i<j
                 }
            }
          if (sum(TRFUM[i,])==0)
            {
             root<-root+1
             root_node[root]<-i
            }
         }
          mismatch<-0      # count the mismatches
          if (length(InOutPair)>0)
          for (i in 1:length(InOutPair))
            {
             if (InOutPair[[i]][length(InOutPair[[i]])]<=10)
               {
                in_node<-c(InOutPair[[i]][2], InOutPair[[i]][3])
                out_node<-InOutPair[[i]][1]
                if (InOutPair[[i]][length(InOutPair[[i]])]==1)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitAnd(GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==2)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], 1-bitAnd(GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==3)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitOr(GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==4)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], 1-bitOr(GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==5)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitOr(1-GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==6)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitOr(GeneData[in_node[1],1:(SampleSize-1)], 1-GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==7)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitAnd(1-GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==8)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitAnd(GeneData[in_node[1],1:(SampleSize-1)], 1-GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==9)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitXor(GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==10)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], 1-bitXor(GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
               }
             if (InOutPair[[i]][length(InOutPair[[i]])]>10)
               {
                in_node<-InOutPair[[i]][2]
                out_node<-InOutPair[[i]][1]
                if (InOutPair[[i]][length(InOutPair[[i]])]==11)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], GeneData[in_node,1:(SampleSize-1)]))
                if (InOutPair[[i]][length(InOutPair[[i]])]==12)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], 1-GeneData[in_node,1:(SampleSize-1)]))
               }
            }
         pseudo_count<-0.0001
         mismatch<-mismatch+pseudo_count   # mismatch may be 0
         Perror<-mismatch/(ii*SampleSize+pseudo_count); #print("Perror="); print(Perror)
         ErrorFactor<-mismatch*log(Perror)+(ii*SampleSize-mismatch+pseudo_count)*log(1-Perror)
         ErrorPrior<-(prior_para[num.node+1,1]-1)*log(Perror)+(prior_para[num.node+1,2]-1)*log(1-Perror)
         RootFactor<-numeric(); RootPrior<-numeric()
         succ_count<-numeric(); succ_prob<-numeric()
         for (i in 1:length(root_node))
           {
            succ_count[i]<-sum(GeneData[root_node[i],])+pseudo_count  # succ_count may be 0
            succ_prob[i]<-succ_count[i]/(SampleSize+pseudo_count)
            RootFactor[i]<-succ_count[i]*log(succ_prob[i])+(SampleSize+pseudo_count-succ_count[i])*log(1-succ_prob[i])

            RootPrior[i]<-(prior_para[root_node[i],1]-1)*log(succ_prob[i])+(prior_para[root_node[i],2]-1)*log(1-succ_prob[i])
           }
         if (length(InOutPair)==0)        # all nodes are root nodes
           {
            likelihood<-sum(RootFactor)
            post_para<-sum(RootFactor) + sum(RootPrior)
            log_post_model<-0
            for (i in 1:nrow(TRFUM))
              {
               nume<-lbeta(prior_para[i,1]+sum(GeneData[i,]), prior_para[i,2]+SampleSize-sum(GeneData[i,])) # lbeta=log(beta)
               deno<-lbeta(prior_para[i,1], prior_para[i,2])
               log_post_model<-log_post_model+nume-deno
              }
            log_post_model=log_post_model+length(TRFUM[TRFUM>0])*log(penalty)
            ErrorFactor<-NA; Perror<-NA; mismatch<-NA
           }
         if (length(InOutPair)>0)        # exist non root nodes
           {
            likelihood<-ErrorFactor+sum(RootFactor)      # this is the likelihood
            post_para<-ErrorFactor+sum(RootFactor) + sum(RootPrior)+ErrorPrior    # this is the posterior of (T, F, theta)
            log_post_model<-0                             # this is the posterior of (T, F)
            for (i in 1:length(root_node))
              {
               index<-root_node[i]
               nume<-lbeta(prior_para[index,1]+sum(GeneData[index,]), prior_para[index,2]+SampleSize-sum(GeneData[index,]))
               deno<-lbeta(prior_para[index,1], prior_para[index,2])
               log_post_model<-log_post_model+nume-deno
              }
            noise_nume<-lbeta(mismatch+prior_para[num.node+1,1], length(InOutPair)*SampleSize-mismatch+prior_para[num.node+1,2])
            noise_deno<-lbeta(prior_para[num.node+1,1], prior_para[num.node+1,2])
            log_post_model<-log_post_model+noise_nume-noise_deno
            log_post_model=log_post_model+length(TRFUM[TRFUM>0])*log(penalty)    # add penalty to the posterior of (T, F)
           }
         result<-list()
         result[[1]]<-c(ErrorFactor, RootFactor, likelihood, post_para, log_post_model)
         para_sample<-rep(NA,num.node)
         for (i in 1:num.node)
          if (i %in% root_node==T)
            for (j in 1:length(root_node))
              if (i==root_node[j])
                para_sample[i]<-succ_prob[j]
         result[[2]]<-c(para_sample, mismatch, Perror)
     return(result)
   }
#############################################################################################################################
#############################################################################################################################
GenerateNetwork<-function(num.node)
{
 loop=1
 while (loop!=0)
  {
   tent_incid_matrix<-matrix(0,nrow=num.node, ncol=num.node)   # define random incidence matrix, ancestor matrix, transition matrix
   tent_trans_matrix<-matrix(0,nrow=num.node, ncol=num.node)
   for (i in 1:num.node)
     {
      u1=runif(1); uu=4       # uu determines the network complexity
      if (u1>1/10 & u1<=uu/10)
       {
        position=sample(seq(num.node)[-i],1); tent_incid_matrix[i,position]=1
        u2=runif(1)
        if (u2>0.5)
          tent_trans_matrix[i,position]=11
        if (u2<0.5)
          tent_trans_matrix[i,position]=12
       }
      if (u1>uu/10)
       {
        position=sample(seq(num.node)[-i],2); tent_incid_matrix[i,position[1]]=1; tent_incid_matrix[i,position[2]]=1
        func=sample(10,1); tent_trans_matrix[i,position[1]]=func; tent_trans_matrix[i,position[2]]=func
       }
     }
   tent_ances_matrix=update.ancestor_matrix(tent_incid_matrix)
   loop=check.ances.matrix(tent_ances_matrix)
  }
  return(tent_trans_matrix)
}
#############################################################################################################################
#############################################################################################################################
GenerateSample<-function(trans_matrix)
{
 node_ances=matrix(nrow=num.node, ncol=2)
 GeneData=matrix(0, nrow=num.node, ncol=SampleSize)
 incid_matrix<-trans_matrix
    for (i in 1:nrow(trans_matrix))
     for (j in 1:ncol(trans_matrix))
      if (trans_matrix[i,j]>0)
       incid_matrix[i,j]<-1
 ances_matrix=update.ancestor_matrix(incid_matrix)
 for (i in 1:num.node)
  {
   node_ances[i,1]=i
   node_ances[i,2]=sum(ances_matrix[i,])
  }
 node_ances=node_ances[order(node_ances[,2]),]
 for (i in 1:nrow(node_ances))
   {
    if (node_ances[i,2]==0)
      GeneData[node_ances[i,1],]=rbinom(SampleSize,1,prob=para[node_ances[i,1]])
    if (node_ances[i,2]!=0)
      {
       parent<-numeric(); ii=1
       for (j in 1:ncol(incid_matrix))
         if (incid_matrix[node_ances[i,1], j]!=0)
           {
            parent[ii]=j; ii=ii+1
           }
       func=trans_matrix[node_ances[i,1],parent[1]]
       if (func==1)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitAnd(GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==2)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(1-bitAnd(GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==3)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitOr(GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==4)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(1-bitOr(GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==5)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitOr(1-GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==6)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitOr(GeneData[parent[1],1:(SampleSize-1)], 1-GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==7)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitAnd(1-GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==8)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitAnd(GeneData[parent[1],1:(SampleSize-1)], 1-GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==9)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitXor(GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==10)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(1-bitXor(GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==11)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(GeneData[parent[1],1:(SampleSize-1)], error[node_ances[i,1],1:(SampleSize-1)])
       if (func==12)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(1-GeneData[parent[1],1:(SampleSize-1)], error[node_ances[i,1],1:(SampleSize-1)])
      }
   }
 return(GeneData)
}
#############################################################################################################################
#############################################################################################################################
BF1<-function(test.stat)    # model g_k=g_i and g_j
  { #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[4],test.stat[6],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
  }
BF2<-function(test.stat)    # model g_k=complement(g_i and g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[3],test.stat[5],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
    return(BIC.value)
 }
BF3<-function(test.stat)    # model g_k=(g_i or g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[3],test.stat[5],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF4<-function(test.stat)    # model g_k=complement(g_i or g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[4],test.stat[6],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF5<-function(test.stat)    # model g_k=complement(g_i) or g_j
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[3],test.stat[6],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF6<-function(test.stat)    # model g_k=g_i or complement( g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[4],test.stat[5],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF7<-function(test.stat)    # model g_k=complement(g_i) and g_j
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[3],test.stat[6],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF8<-function(test.stat)    # model g_k=g_i and complement(g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[4],test.stat[5],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF9<-function(test.stat)    # model g_k=g_i xor g_j
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[3],test.stat[5],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
    return(BIC.value)
 }
 BF10<-function(test.stat)    # model g_k=complement(g_i xor g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[4],test.stat[6],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF11<-function(test.stat)    # model g_k=g_i
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[4],test.stat[5],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF12<-function(test.stat)    # model g_k=g_j
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[3],test.stat[6],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF13<-function(test.stat)    # model g_k=complement(g_i)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[3],test.stat[6],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
 BF14<-function(test.stat)    # model g_k=complement( g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[4],test.stat[5],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
 ProposalConstruction<-function(GeneData)
 {
gene.data=GeneData
sample.size<-ncol(gene.data)
num.node<-nrow(gene.data)
error.prop<-0.4; pseudo.count<-0.01
sample.size<-sample.size+pseudo.count*8; threshold<-sample.size*error.prop
candidate.prior<-list()
kk<-1
for (i in 1: nrow(gene.data))
  for (j in 1: nrow(gene.data))
   if (j!=i)
   for (k in 1: nrow(gene.data))
    if (k!=i & k!=j)
      {  # print(i)
       test.result<-list()  # save test results for all possible relations
       gene.triplet<-rbind(rbind(gene.data[i,1:(SampleSize-1)], gene.data[j,1:(SampleSize-1)]),gene.data[k,2:SampleSize])
        c000<-0; c001<-0; c010<-0; c011<-0; c100<-0; c101<-0; c110<-0; c111<-0
        for (ii in 1: (ncol(gene.triplet)-1))   # counts in each cell
     {
      if (gene.triplet[1,ii]==0 & gene.triplet[2,ii]==0 & gene.triplet[3,(ii+1)]==0)
         c000<-c000+1
      if (gene.triplet[1,ii]==0 & gene.triplet[2,ii]==0 & gene.triplet[3,(ii+1)]==1)
         c001<-c001+1
      if (gene.triplet[1,ii]==0 & gene.triplet[2,ii]==1 & gene.triplet[3,(ii+1)]==0)
         c010<-c010+1
      if (gene.triplet[1,ii]==0 & gene.triplet[2,ii]==1 & gene.triplet[3,(ii+1)]==1)
         c011<-c011+1
      if (gene.triplet[1,ii]==1 & gene.triplet[2,ii]==0 & gene.triplet[3,(ii+1)]==0)
         c100<-c100+1
      if (gene.triplet[1,ii]==1 & gene.triplet[2,ii]==0 & gene.triplet[3,(ii+1)]==1)
         c101<-c101+1
      if (gene.triplet[1,ii]==1 & gene.triplet[2,ii]==1 & gene.triplet[3,(ii+1)]==0)
         c110<-c110+1
      if (gene.triplet[1,ii]==1 & gene.triplet[2,ii]==1 & gene.triplet[3,(ii+1)]==1)
         c111<-c111+1
      }
       test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)   #  generate random sample
       test.result[[1]]<-c(i,j, k, 1, BF1(test.stat))
       test.result[[2]]<-c(i,j, k, 2, BF2(test.stat))
       test.result[[3]]<-c(i,j, k, 3, BF3(test.stat))
       test.result[[4]]<-c(i,j, k, 4, BF4(test.stat))
       test.result[[5]]<-c(i,j, k, 5, BF5(test.stat))
       test.result[[6]]<-c(i,j, k, 6, BF6(test.stat))
       test.result[[7]]<-c(i,j, k, 7, BF7(test.stat))
       test.result[[8]]<-c(i,j, k, 8, BF8(test.stat))
       test.result[[9]]<-c(i,j, k, 9, BF9(test.stat))
       test.result[[10]]<-c(i,j, k, 10, BF10(test.stat))
       test.result[[11]]<-c(i,j, k, 11, BF11(test.stat)) # model g_k=g_i
       test.result[[12]]<-c(i,j, k, 12, BF12(test.stat)) # model g_k=g_j
       test.result[[13]]<-c(i,j, k, 13, BF13(test.stat)) # model g_k=complement(g_i)
       test.result[[14]]<-c(i,j, k, 14, BF14(test.stat))  # model g_k=complement(g_j)
       # save the most likely of all 14 relations by their false counts.
       # *******************here it may filter out some true ones due to noise and model uncertainty ****************
       miscount<-numeric(); jj<-1
       for (ii in 1:length(test.result))
         if (length(test.result[[ii]])==5)
          {
            miscount[jj]<-test.result[[ii]][5]
            jj<-jj+1
          }

         if (length(miscount)>0)
         {
            min.miscount<-min(miscount)
           for (ii in 1:length(test.result))
           if ( length(test.result[[ii]])==5 & test.result[[ii]][5]==min.miscount)
             {
              candidate.prior[[kk]]<-test.result[[ii]]
              kk<-kk+1
             }
         }   # for each pariwise genes or gene triplet, only the most likely one is output
    }
candidate<-matrix(nrow=length(candidate.prior), ncol=length(candidate.prior[[1]]))
for ( i in 1:length(candidate.prior))
  candidate[i,]<-candidate.prior[[i]]
order.candidate<-candidate[order(candidate[,3]),]   #order by output variables
##################################   # transform BIC into proportion
 num.node<-max(order.candidate[,3])
 triplet<-list();j1<-1; pairwise<-list(); j2<-1
  for (j in 1:nrow(order.candidate))
   {
   if (order.candidate[j,4]<=10)
     {
       triplet[[j1]]<-order.candidate[j,]
       j1<-j1+1
     }
   if (order.candidate[j,4]>10)
     {
      if (order.candidate[j,4]==11)
      {
        pairwise[[j2]]<-c(order.candidate[j,1], order.candidate[j,3], 11, order.candidate[j,5])
        j2<-j2+1
      }
       if (order.candidate[j,4]==12)
      {
        pairwise[[j2]]<-c(order.candidate[j,2], order.candidate[j,3], 11, order.candidate[j,5])  # model 11, g_k=g_i
        j2<-j2+1
      }
       if (order.candidate[j,4]==13)
      {
        pairwise[[j2]]<-c(order.candidate[j,1], order.candidate[j,3], 12, order.candidate[j,5])
        j2<-j2+1
      }
      if (order.candidate[j,4]==14)
      {
        pairwise[[j2]]<-c(order.candidate[j,2], order.candidate[j,3], 12, order.candidate[j,5])# model 12: complement relation
        j2<-j2+1
      }
     }
   }
candidate.triplet<-matrix(nrow=length(triplet),ncol=5); weighted.triplet<-matrix(nrow=length(triplet),ncol=6)
candidate.pairwise<-matrix(nrow=length(pairwise),ncol=4)
for (i in 1:length(triplet))
  {
   candidate.triplet[i,]<-triplet[[i]]
   for (j in 1:5)
    weighted.triplet[i,j]=candidate.triplet[i,j]
   weighted.triplet[i,6]=candidate.triplet[i,5]
  }
for (i in 1:length(pairwise))
  candidate.pairwise[i,]<-pairwise[[i]]
 unique.pairwise<-unique(candidate.pairwise)
 weighted.pairwise<-matrix(nrow=length(unique.pairwise),ncol=5)
 weighted.pairwise=cbind(unique.pairwise,unique.pairwise[,4])
constant<-0.001
for (i in 1:num.node)
 if (i %in% candidate.triplet[,3]==T)
 for (j in 1:nrow(candidate.triplet))
   if (candidate.triplet[j,3]==i)
       {
        trip.matrix<-matrix(ncol=5)
        trip.matrix<-candidate.triplet[candidate.triplet[,3]==i,1:5]
        trip.matrix<-data.matrix(trip.matrix)  # avoid pair.matrix to be "numeric"
        if (ncol(trip.matrix)==1)
         trip.matrix<-t(trip.matrix)
        score<-1/trip.matrix[,5]          # use reciprocal of miscunt as weight
        prop<-(score+constant)/(sum(score+constant))
        trip.matrix[,5]<-prop
        weighted.triplet[weighted.triplet[,3]==i,]<-cbind(trip.matrix, weighted.triplet[weighted.triplet[,3]==i,6])
        }

for (i in 1:num.node)
  if (i %in% unique.pairwise[,2]==T)
  for (j in 1:nrow(unique.pairwise))
   if (unique.pairwise[j,2]==i)
       {
        pair.matrix<-matrix(ncol=4)
        pair.matrix<-unique.pairwise[unique.pairwise[,2]==i,1:4]
        pair.matrix<-data.matrix(pair.matrix)  # avoid pair.matrix to be "numeric"
        if (ncol(pair.matrix)==1)
        pair.matrix<-t(pair.matrix)
        score<-1/(pair.matrix[,4]+constant)      # note here
        prop<-score/(sum(score))
        pair.matrix[,4]<-prop
        weighted.pairwise[weighted.pairwise[,2]==i,]<-cbind(pair.matrix, weighted.pairwise[weighted.pairwise[,2]==i,5])
        }
CandidateTriplet=weighted.triplet
CandidatePairwise=weighted.pairwise
Candidate=list()
Candidate[[1]]=CandidateTriplet
Candidate[[2]]=CandidatePairwise
return(Candidate)
 }
 ###################################################################################################
ConstructInitial<-function(Candidate)
{
prior.triple=Candidate[[1]]
prior.pairwise=Candidate[[2]]
trans_func_matrix<-matrix(0,nrow=num.node,ncol=num.node)
prop_trans_func_matrix<-matrix(0,nrow=num.node, ncol=num.node)
prop_incid_matrix<-matrix(0, nrow=num.node, ncol=num.node)
prop_ances_matrix<-matrix(0, nrow=num.node, ncol=num.node)
incid_matrix<-matrix(0,nrow=num.node,ncol=num.node)
ratio<-0.5
for (i in 1:num.node)
  {
      prop_trans_func_matrix<-trans_func_matrix
      prop_incid_matrix<-incid_matrix
   #################################
      pairwise.prior.set<-matrix()        # clear variables     # determine candidate parents for  node i
      pairwise.prior.pare<-matrix()
      pairwise.prior.set<-prior.pairwise[prior.pairwise[,2]==i,]    # priors for node update_order[k]
      pairwise.prior.set<-data.matrix(pairwise.prior.set)
      if (ncol(pairwise.prior.set)==1)
       {
         pairwise.prior.set<-t(pairwise.prior.set)
         pairwise.prior.pare<-t(data.matrix(pairwise.prior.set[,1]))
       }
      if (ncol(pairwise.prior.set)>1 & nrow(pairwise.prior.set)>0)
        pairwise.prior.pare<-pairwise.prior.set[,1]

      triplet.prior.set<-matrix()
      triplet.prior.pare<-matrix()
      triplet.prior.set<-prior.triplet[prior.triplet[,3]==i,]
      triplet.prior.set<-data.matrix(triplet.prior.set)
      if (ncol(triplet.prior.set)>1 & nrow(triplet.prior.set)>0)
        triplet.prior.pare<-triplet.prior.set[,1:2]
      if (ncol(triplet.prior.set)==1)
        {
          triplet.prior.set<-t(triplet.prior.set)
          triplet.prior.pare<-t(data.matrix(triplet.prior.set[,1:2]))
         }
      prop_prob<-runif(1); aa<-0
      ##############################################
      if (prop_prob>=ratio)   # consider two parents
       {
        if (nrow(triplet.prior.set)>=1)
                 {
                    candidate.pare.set<-list()
                    for (jj in 1:nrow(triplet.prior.set))
                      if (length(triplet.prior.set[jj,1:2])==2)
                        {
                         aa<-aa+1
                         candidate.pare.set[[aa]]<-triplet.prior.set[jj,]
                        }
                    if (length(candidate.pare.set)>0)
                        {
                         score<-numeric()
                         for (jj in 1:length(candidate.pare.set))
                           score[jj]<-candidate.pare.set[[jj]][5]
                         bbb<-matrix(nrow=length(candidate.pare.set),ncol=6)
                         for (jjj in 1:length(candidate.pare.set))
                           bbb[jjj,]<-candidate.pare.set[[jjj]]

                         max.score.candidate<-bbb[bbb[,5]==max(score),]
                         if (is.numeric(max.score.candidate)==T)
                           add_parent=max.score.candidate
                         if (is.matrix(max.score.candidate)==T)  
                         if (nrow(max.score.candidate)>1)
                           add_parent<-max.score.candidate[sample.int(nrow(max.score.candidate),1),]

                         add_two_parent<-c(add_parent[1],add_parent[2])
                      }
                 }
           if (aa>0)
             {
              func_order<-add_parent[4]
              prop_trans_func_matrix[i, add_two_parent[1]]<-func_order
              prop_trans_func_matrix[i, add_two_parent[2]]<-func_order
              for (ii in 1:nrow(prop_trans_func_matrix))
                for (jj in 1:ncol(prop_trans_func_matrix))
                 if (prop_trans_func_matrix[ii,jj]>0)
                   prop_incid_matrix[ii,jj]<-1
              prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
              if (check.ances.matrix(prop_ances_matrix)==0)
               {
                trans_func_matrix[i,]<-prop_trans_func_matrix[i,]
                incid_matrix[i,]<-prop_incid_matrix[i,]
               }
             }
          }
      ##########################################
      if (prop_prob<ratio) # consider one parent
       if (nrow(pairwise.prior.set)>0)
         {
            candidate.pare.set<-pairwise.prior.set
            if (is.numeric(candidate.pare.set)==T)
              add_parent<-candidate.pare.set
            if (is.matrix(candidate.pare.set)==T)
              add_parent<-candidate.pare.set[candidate.pare.set[,4]==max(candidate.pare.set[,4]),]
            if (is.matrix(add_parent)==T)
              if ( nrow(add_parent)>1)
                add_parent<-add_parent[sample.int(nrow(add_parent),1),]

            add_one_parent<-add_parent[1]

           func_order<-add_parent[3]
           prop_trans_func_matrix[i, add_one_parent]<-func_order
           for (ii in 1:nrow(prop_trans_func_matrix))
             for (jj in 1:ncol(prop_trans_func_matrix))
               if (prop_trans_func_matrix[ii,jj]>0)
                 prop_incid_matrix[ii,jj]<-1
           prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
           if (check.ances.matrix(prop_ances_matrix)==0)
            {
             trans_func_matrix[i,]<-prop_trans_func_matrix[i,]
             incid_matrix[i,]<-prop_incid_matrix[i,]
            }
          }
  }
  return(trans_func_matrix)
}
#######################################################################################
###############  MCMC
prop.ratio<-0.1       # proposal information is used with probability prop.ratio
prop_beta1<-2; prop_beta2<-100
Trans_Func<-seq(1,12)     # all transition functions where each node has two parents at most
num.node=20; SampleSize=200
error.prop<-0.2; pseudo.count<-0.01
threshold<-SampleSize*error.prop
num_update=200
penalty=0.1
 prior_para<-matrix(3, nrow=(num.node+1), ncol=2)
 prior_para[num.node+1,1]<-2;prior_para[num.node+1,2]<-100    # for error, it should be beta1e<beta2e

 para<-numeric(); ncp=0
 for (i in 1:nrow(prior_para))
    para[i]<-rbeta(1, prior_para[i,1], prior_para[i,2],ncp)
 para[21]=0.1
 error<-matrix(0, nrow=num.node,ncol=SampleSize)
 for(i in 1:num.node)
  error[i,]<-rbinom(SampleSize,1,prob=para[num.node+1])
###############
  true_network=GenerateNetwork(num.node)       # randomly generate a network 

  GeneData=matrix(nrow=num.node, ncol=SampleSize)   # based on the generated network, create a data set 
  GeneData=GenerateSample(true_network)
 # GeneData<-read.table("GeneData_randomdatanoise0.1proposalinitialpenalty0.1multiplechains_3.txt",header=T)  #true network
 # GeneData<-data.matrix(GeneData)
 # true_network<-read.table("TrueNetwork_randomdatanoise0.1proposalinitialpenalty0.1multiplechains_3.txt",header=T)
 # true_network<-data.matrix(true_network)
  true_incid_matrix<-true_network
    for (i in 1:nrow(true_incid_matrix))
     for (j in 1:ncol(true_incid_matrix))
      if (true_incid_matrix[i,j]>0)
       true_incid_matrix[i,j]<-1
  true_logpost=Error_LLH(true_network)[[1]][length(Error_LLH(true_network))]

  Candidate=ProposalConstruction(GeneData)    # create the proposal for generated data  
  prior.triplet<-Candidate[[1]]
 # prior.triplet=read.table("proposal_triplet_true_network.txt",header=T)
  prior.pairwise<-Candidate[[2]]
 # prior.pairwise=read.table("proposal_pairwise_true_network.txt",header=T)
  #####################################
  StandAlone=matrix(nrow=num.node,ncol=2)   # calculate the BIC for each node
   for (i in 1:num.node)
     {
      succ.prob=(sum(GeneData[i,])+pseudo.count)/(SampleSize+pseudo.count)  # calculate BIC for standalone gene 
      BIC.value=-2*(sum(GeneData[i,]+pseudo.count)*log(succ.prob)+(SampleSize-sum(GeneData[i,]))*log(1-succ.prob))+1*log(SampleSize)
      post.data=exp(-0.5*BIC.value)
      StandAlone[i,1]=i; StandAlone[i,2]=BIC.value
     }
  #StandAlone=read.table("proposal_standalone_true_network.txt",header=T)
  ###################################  multiple independent chains
#  All_Trans_Func_Matrix=list()
#  All_Correct_Rate=list()
#  All_Logpost=list()
  ###################################

 for (iii in 1:1)       # iii: number of simulations 
  {
  trans_func_matrix=ConstructInitial(Candidate)    # use the randomly selected initial 
  #trans_func_matrix=read.table("correctrate0.4trans_func_matrix_twostepneighborhood3_8_works_finalnetwork.txt",header=T) # use specific initial 
  # trans_func_matrix=read.table("CorrectRate0.55AsStartingNetwork.txt",header=T)
  incid_matrix<-trans_func_matrix
    for (i in 1:nrow(incid_matrix))
     for (j in 1:ncol(incid_matrix))
      if (incid_matrix[i,j]>0)
       incid_matrix[i,j]<-1
  ances_matrix=update.ancestor_matrix(incid_matrix)


  Incidence_Matrix<-list()
  Ancestor_Matrix<-list()     # Matrix: Ancestor Matrix recording ancesor-offspring relatons for the whole chain
  Trans_Func_Matrix<-list()    # Matrix: Transition Function Matrix for the whole chain
  Sample_Matrix<-matrix(nrow=num.node*num_update+1, ncol=num.node+2)
  Sample_Matrix[1,]<-Error_LLH(trans_func_matrix)[[2]]
  Incidence_Matrix[[1]]<-incid_matrix
  Ancestor_Matrix[[1]]<-ances_matrix
  Trans_Func_Matrix[[1]]<-trans_func_matrix
  num<-numeric();logpost<-numeric()
  all_logpost<-numeric(); aa<-Error_LLH(Trans_Func_Matrix[[1]]);  all_logpost[1]<-aa[[1]][length(aa[[1]])]
  logpost[1]<-all_logpost[1]
  n<-1; num[1]<-1;  iter<-1 ; jump_point<-numeric(); jump_point[1]<-1   # paramters for each chain
for (ii in 1: num_update)       # run ii full rounds, with each round of num.node times
{
 if (ii<=round(0.1*num_update))       # use adaptive ratio of using proposal information 
   prop.ratio=0.1 
 if (ii>round(0.1*num_update))
   prop.ratio=0.9 
     
 update_order<-sample.int(num.node,num.node,replace=FALSE)
 for (k in 1: length(update_order))         #   consider the updating node g_k
  {
    cat("ii", ii,  "th iteration is running", "\n")
    cat("n=", n,"\n" )
    old<-n
 #   print("iter="); print(iter)
    cat("true_logpost=", true_logpost, "\n")
    cat("all_logpost=", all_logpost[iter],"\n")
    current_incid_matrix<-Incidence_Matrix[[iter]]
    current_ances_matrix<-Ancestor_Matrix[[iter]]
    current_trans_func_matrix<-Trans_Func_Matrix[[iter]]
    current_post<-all_logpost[iter]   # here current_post is a scale

    parent_of_update<-numeric();j<-1      # find the parent for node update_order[k]
    for (i in 1: ncol(current_incid_matrix))
      if (current_incid_matrix[update_order[k],i]!=0)
        {
         parent_of_update[j]<-i; j<-j+1
        }
    swap_candi<-numeric(); j2<-1       # swap_candi  is used for swapping parent action     # swap_candi & legal_parent depend on current_parent
    legal_parent<-numeric(); j1<-1    # legal_parent is used for adding parent action
    for (i in 1:num.node)
     if ( i!=update_order[k] & current_ances_matrix[i,update_order[k]]!=1 & (i  %in% parent_of_update == FALSE))
      {
       legal_parent[j1]<-i; j1<-j1+1 #; print(i)
      }
      swap_candi<-legal_parent
     #################################    # determine the proposal for current updating node
      pairwise.prior.set<-matrix()        # clear variables
      pairwise.prior.pare<-matrix()
      pairwise.prior.set<-prior.pairwise[prior.pairwise[,2]==update_order[k],]    # pairwise proposal for node update_order[k]
      pairwise.prior.set<-data.matrix(pairwise.prior.set)
      if (ncol(pairwise.prior.set)==1)
       {
         pairwise.prior.set<-t(pairwise.prior.set)
         pairwise.prior.pare<-t(data.matrix(pairwise.prior.set[,1]))
    #     all.maxscore.pairwise=pairwise.prior.set
       }
      if (ncol(pairwise.prior.set)>1 & nrow(pairwise.prior.set)>0)
       {
         pairwise.prior.pare<-pairwise.prior.set[,1]
    #     all.maxscore.pairwise=pairwise.prior.set[pairwise.prior.set[,4]==max(pairwise.prior.set[,4]),]
       }  

      triplet.prior.set<-matrix()
      triplet.prior.pare<-matrix()
      triplet.prior.set<-prior.triplet[prior.triplet[,3]==update_order[k],]      # triplet proposal for node update_order[k]
      triplet.prior.set<-data.matrix(triplet.prior.set)
      if (ncol(triplet.prior.set)>1 & nrow(triplet.prior.set)>0)
        {
         triplet.prior.pare<-triplet.prior.set[,1:2]
   #      all.maxscore.triplet=triplet.prior.set[triplet.prior.set[,5]==max(triplet.prior.set[,5]),]
        } 
      if (ncol(triplet.prior.set)==1)
        {
          triplet.prior.set<-t(triplet.prior.set)
          triplet.prior.pare<-t(data.matrix(triplet.prior.set[,1:2]))
   #       all.maxscore.triplet=triplet.prior.set
         }
      ###########################  find the minimum mismatches for gene pairs and gene triplet
   #   BIC.pairwise=0; BIC.triplet=0 # case1:  "numeric case"
   #     if(is.numeric(all.maxscore.pairwise)==T & length(all.maxscore.pairwise)>0)
   #        maxscore.pairwise=all.maxscore.pairwise
   #       
   #     if(is.numeric(all.maxscore.triplet)==T & length(all.maxscore.triplet)>0)
   #        maxscore.triplet=all.maxscore.triplet
   #       
   #     if (is.matrix(all.maxscore.pairwise)==T)
   #      if (nrow(all.maxscore.pairwise)>0)             # randomly select one row if maxscore.pairwsie has multiple rows
   #       if (is.matrix(all.maxscore.triplet)==T)
   #        if(nrow(all.maxscore.triplet)>0)
   #      {
   #       maxscore.pairwise=maxscore.pairwise[sample.int(nrow(maxscore.pairwise),1),]
   #       maxscore.triplet=maxscore.triplet[sample.int(nrow(maxscore.triplet),1),]
   #      }
   #      BIC.pairwise=maxscore.pairwise[5]
   #      BIC.triplet=maxscore.triplet[6]
   #      BIC.standalone=StandAlone[update_order[k],2]
    # pairwise.ratio=BIC.pairwise/(BIC.pairwise+BIC.triplet)
    pairwise.ratio=1/2
############&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# case one: no parent
      if (length(parent_of_update)==0) #  there are only add-parent moves, every time only one move can be proposed
        {
         num_legal_parent<-length(legal_parent)
         if (num_legal_parent>1)          
            uu<-runif(1)
         if (num_legal_parent==0)       # if no candiate nodes are available for adding parent, then ignore this node
            uu<-0
         if (num_legal_parent==1||(num_legal_parent>1 & uu>=pairwise.ratio))
           {
             ################################     Every move should have a reasonable proposal probability
             # proposal move 1: add one parent
            prop.legal.overlap<-intersect(pairwise.prior.pare, legal_parent)

            prop.prob<-runif(1); aa<-0
            if (prop.prob>=prop.ratio)    # proposal information is used with probability 1-prop.ratio
             if (length(prop.legal.overlap)>=1 & nrow(pairwise.prior.set)>0)
               {
                 aa<-aa+1
                 candidate.pare.set<-pairwise.prior.set[pairwise.prior.set[,1]%in%prop.legal.overlap==T,]
                 if (is.numeric(candidate.pare.set)==T)
                    add_parent<-candidate.pare.set
                 if (is.matrix(candidate.pare.set)==T)        # should modify to use multinomial distribution
                    add_parent<-candidate.pare.set[candidate.pare.set[,4]==max(candidate.pare.set[,4]),]  #use the most likely one
                 if (is.matrix(add_parent)==T)
                  if ( nrow(add_parent)>1)
                    add_parent<-add_parent[sample.int(nrow(add_parent),1),]

                 add_one_parent<-add_parent[1]
               }
             if (prop.prob<prop.ratio||aa==0)
              {
                 sample_node<-sample.int(length(legal_parent),1,replace=F)
                 add_one_parent<-legal_parent[sample_node]
              }
             prop_incid_matrix<-current_incid_matrix; prop_incid_matrix[update_order[k], add_one_parent]<-1
             prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
             prop_trans_func_matrix<-current_trans_func_matrix
             if (aa>0)
               prop_trans_func_matrix[update_order[k], add_one_parent]<-add_parent[3]     # function prior
             if (aa==0)
               {
                func_order<-10+sample.int(2,1)
                prop_trans_func_matrix[update_order[k], add_one_parent]<-func_order
               }
             add_one_prob<-1/length(legal_parent)
             prop_trans_func_prob<-1/2  #  there are in all two types of boolean functions to choose for pairwise genes
             xxx<-Error_LLH(prop_trans_func_matrix)
             prop_sample<-xxx[[2]]
             prop_post<-xxx[[1]][length(xxx[[1]])]
             prop_sample_prob<-numeric(); p2c_prob<-numeric()   # calculate acceptance probability
             curr_sample_prob<-numeric(); c2p_prob<-numeric()
             prop_sample_prob[1]<-prop_post;prop_sample_prob[2]<-0      #  this is the prior & likelihood
             curr_sample_prob[1]<-current_post[1];curr_sample_prob[2]<-0
             prop_sample_prob[3]<-log(1/2); curr_sample_prob[3]<-log(1)    #  this is P(R|T)

             if (num_legal_parent>1) # Q(T_c|T_p)
              p2c_prob[1]<-1/3       # there are 3 possible moves, add one; remove one; swap one
             if (num_legal_parent==1)
              p2c_prob[1]<-1

             if (num_legal_parent>1)
                 c2p_prob[1]<-1/2*add_one_prob   # this is Q(T_p|T_c) the probability of move type and specific node involved 
             else
                 c2p_prob[1]<-add_one_prob
             p2c_prob[2]<-1   # Q(R_c|T_c)
             c2p_prob[2]<-prop_trans_func_prob  # Q(R_p|T_p)
             p2c_prob[3]<-1 # Q(\theta_c|T_c, R_c), i.e. the prior probability for new root node
             c2p_prob[3]<-1                     # Q(\theta_p|T_p, R_p)
                                                                                  
             nume<-sum((prop_sample_prob))+sum(log(p2c_prob))
             deno<-sum((curr_sample_prob))+sum(log(c2p_prob))
             acce_prob<-exp(nume-deno)
             ratio<-runif(1)
             if (ratio<=acce_prob)
                 n<-n+1
           }
            ##################################
            # proposal move 2 add two parents
            if(num_legal_parent>1 & uu<pairwise.ratio)      # add two parents one time
             {
              prop.prob<-runif(1); aa<-0
              if(prop.prob>=prop.ratio)
                if (nrow(triplet.prior.set)>=1)
                 {
                    candidate.pare.set<-list()
                    for (jj in 1:nrow(triplet.prior.set))
                      if (length(intersect(triplet.prior.set[jj,1:2], legal_parent))==2)
                        {
                         aa<-aa+1
                         candidate.pare.set[[aa]]<-triplet.prior.set[jj,]
                        }
                    if (length(candidate.pare.set)>0)
                        {
                         score<-numeric()
                         for (jj in 1:length(candidate.pare.set))
                           score[jj]<-candidate.pare.set[[jj]][5]
                         bbb<-matrix(nrow=length(candidate.pare.set),ncol=6)
                         for (jjj in 1:length(candidate.pare.set))
                           bbb[jjj,]<-candidate.pare.set[[jjj]]

                         max.score.candidate<-bbb[bbb[,5]==max(score),]
                         if (is.numeric(max.score.candidate)==T)
                           add_parent<-max.score.candidate
                         if (is.matrix(max.score.candidate)==T)
                          if (nrow(max.score.candidate)>1)
                           add_parent<-max.score.candidate[sample.int(nrow(max.score.candidate),1),]

                         add_two_parent<-c(add_parent[1],add_parent[2])
                      }
                   }
               if (prop.prob<prop.ratio||aa==0)
               {
                   sample_two_parent<-sample.int(length(legal_parent),2,replace=FALSE)
                   add_two_parent<-c(legal_parent[sample_two_parent[1]],legal_parent[sample_two_parent[2]])
               }
               prop_incid_matrix<-current_incid_matrix
               prop_incid_matrix[update_order[k], add_two_parent[1]]<-1;prop_incid_matrix[update_order[k], add_two_parent[2]]<-1
               prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
               prop_trans_func_matrix<-current_trans_func_matrix
               if (aa>0)
                 {
                  prop_trans_func_matrix[update_order[k], add_two_parent[1]]<-add_parent[4]
                  prop_trans_func_matrix[update_order[k], add_two_parent[2]]<-add_parent[4]
                 }
               if (aa==0)
                 {
                  func_order<-sample.int(10,1)
                  prop_trans_func_matrix[update_order[k], add_two_parent[1]]<-func_order
                  prop_trans_func_matrix[update_order[k], add_two_parent[2]]<-func_order
                 }
               sample_two_prob<-1/choose(num_legal_parent,2)
               prop_trans_func_prob<-1/10
               xxx<-Error_LLH(prop_trans_func_matrix)
               prop_sample<-xxx[[2]]
               prop_post<-xxx[[1]][length(xxx[[1]])]
               prop_sample_prob<-numeric(); p2c_prob<-numeric() # calculate acceptance probability
               curr_sample_prob<-numeric(); c2p_prob<-numeric()
               prop_sample_prob[1]<-prop_post;prop_sample_prob[2]<-0       #  this is the prior & likelihood
               curr_sample_prob[1]<-current_post[1];curr_sample_prob[2]<-0
               prop_sample_prob[3]<-log(prop_trans_func_prob); curr_sample_prob[3]<-log(1)    #  this is P(R|T)

               if (length(swap_candi)>1)
                 p2c_prob[1]<-1/4 #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding two parents, depending on the number of nodes to swap
               if (length(swap_candi)==1)
                 p2c_prob[1]<-1/3
               if (length(swap_candi)==0)
                 p2c_prob[1]<-1/2
               c2p_prob[1]<-1/2*sample_two_prob   # this is Q(T_p|T_c) there are adding one parent, adding two parents, two moves
               p2c_prob[2]<-1   # Q(R_c|T_c)
               c2p_prob[2]<-prop_trans_func_prob  # Q(R_p|T_p)
               p2c_prob[3]<-1# Q(\theta_c|T_c, R_c), i.e. the prior probability for new root node
               c2p_prob[3]<-1                     # Q(\theta_p|T_p, R_p)

               nume<-sum((prop_sample_prob))+sum(log(p2c_prob))
               deno<-sum((curr_sample_prob))+sum(log(c2p_prob))
               acce_prob<-exp(nume-deno)
               ratio<-runif(1)
               if (ratio<=acce_prob)
                   n<-n+1
             }
         }
##########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# case 2: one parent
       if (length(parent_of_update)==1)
         {
          uu<-runif(1)
          #one.parent.ratio=2/5*(1/BIC.triplet/(1/BIC.triplet+1/BIC.standalone))
          one.parent.ratio=1/6
          #################################
          #proposal move 1: add one parent
           if (length(legal_parent)>0 & uu<one.parent.ratio)
            {
             prop.prob<-runif(1);aa<-0
             if(prop.prob>=prop.ratio)
               if (nrow(triplet.prior.set)>=1 & ncol(triplet.prior.set)>1)
                 {
                    candidate.pare.set<-list()
                    for (jj in 1:nrow(triplet.prior.set))
                      if (length(intersect(parent_of_update,triplet.prior.pare[jj,] ))==1)
                      if (parent_of_update %in% triplet.prior.pare[jj,]==T & setdiff(triplet.prior.pare[jj,], parent_of_update) %in% legal_parent==T)
                        {
                          aa<-aa+1
                          candidate.pare.set[[aa]]<-triplet.prior.set[jj,]
                         }
                    if (length(candidate.pare.set)>0)
                      {
                       score<-numeric()
                       for (jj in 1:length(candidate.pare.set))
                         score[jj]<-candidate.pare.set[[jj]][5]
                        bbb<-matrix(nrow=length(candidate.pare.set),ncol=6)
                        for (jjj in 1:length(candidate.pare.set))
                         bbb[jjj,]<-candidate.pare.set[[jjj]]

                        max.score.candidate<-bbb[bbb[,5]==max(score),]
                        if (is.numeric(max.score.candidate)==T)
                          add_parent<-max.score.candidate
                        if (is.matrix(max.score.candidate)==T)
                         if (nrow(max.score.candidate)>1)
                         add_parent<-max.score.candidate[sample.int(nrow(max.score.candidate),1),]

                        add_one_parent<-setdiff(add_parent[1:2], parent_of_update)
                     }
                  }
               if (prop.prob<prop.ratio||aa==0)
               {
                  sample_node<-sample.int(length(legal_parent),1,replace=FALSE)
                  add_one_parent<-legal_parent[sample_node]
               }
               prop_incid_matrix<-current_incid_matrix ; prop_incid_matrix[update_order[k], add_one_parent]<-1
               prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
               prop_trans_func_matrix<-current_trans_func_matrix
               if (aa>0)
                 {
                  prop_trans_func_matrix[update_order[k], add_one_parent]<-add_parent[4]
                  prop_trans_func_matrix[update_order[k], parent_of_update]<-add_parent[4]
                 }
               if (aa==0)
                 {
                  func_order<-sample.int(10,1)
                  prop_trans_func_matrix[update_order[k], add_one_parent]<-func_order
                  prop_trans_func_matrix[update_order[k], parent_of_update]<-func_order
                 }
               add_one_prob<-1/length(legal_parent)
               prop_trans_func_prob<-1/10  # there are 10 functions to choose
               xxx<-Error_LLH(prop_trans_func_matrix)
               prop_sample<-xxx[[2]]
               prop_post<-xxx[[1]][length(xxx[[1]])]
               prop_sample_prob<-numeric(); p2c_prob<-numeric() # calculate acceptance probability
               curr_sample_prob<-numeric(); c2p_prob<-numeric()
               prop_sample_prob[1]<-prop_post;prop_sample_prob[2]<-0       #  this is the prior & likelihood
               curr_sample_prob[1]<-current_post[1];curr_sample_prob[2]<-0
               prop_sample_prob[3]<-log(prop_trans_func_prob); curr_sample_prob[3]<-log(1)    #  this is P(R|T)

               if (length(swap_candi)>1)
                 p2c_prob[1]<-1/4 #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
               if (length(swap_candi)==1)
                 p2c_prob[1]<-1/3
               if (length(swap_candi)==0)
                 p2c_prob[1]<-1/2

               if (length(swap_candi)>0)
                c2p_prob[1]<-1/3*add_one_prob   # this is Q(T_p|T_c) probability of which move and which node selected
               if (length(swap_candi)==0)
                c2p_prob[1]<-1/2*add_one_prob

               p2c_prob[2]<-1   # Q(R_c|T_c)
               c2p_prob[2]<-prop_trans_func_prob  # Q(R_p|T_p)
               p2c_prob[3]<-1# Q(\theta_c|T_c, R_c), i.e. the prior probability for new root node
               c2p_prob[3]<-1                     # Q(\theta_p|T_p, R_p)

               nume<-sum((prop_sample_prob))+sum(log(p2c_prob))
               deno<-sum((curr_sample_prob))+sum(log(c2p_prob))
               acce_prob<-exp(nume-deno); ratio<-runif(1)
               if (ratio<=acce_prob)
                  n<-n+1
          }
          ##################################
          # proposal move 2, swap one parent
           if (length(swap_candi)>0 & uu>=1/6 & uu<2/6)
             {
              swap.prior.overlap<-intersect(swap_candi, pairwise.prior.pare)
              prop.prob<-runif(1); aa<-0
              if(prop.prob>prop.ratio)
               if (length(swap.prior.overlap)>0 & nrow(pairwise.prior.set)>0)
                {
                 aa<-aa+1
                 swap.candidate<-pairwise.prior.set[pairwise.prior.set[,1]%in%swap.prior.overlap==T,]
                 if (is.numeric(swap.candidate)==T)   # i.e. swap.candidate has one row
                      swap_parent<-swap.candidate
                 if (is.matrix(swap.candidate)==T)
                     swap_parent<-swap.candidate[swap.candidate[,4]==max(swap.candidate[,4]),] # multiple norw may have same weight

                 if (is.matrix(swap_parent)==T)
                 if ( nrow(swap_parent)>1)
                   swap_parent<-swap_parent[sample.int(nrow(swap_parent),1),]  # every candidate has equal chance
                 swap_one_node<-swap_parent[1]
                }
               if (prop.prob<prop.ratio||aa==0)
               {
                 sample_node<-sample.int(length(swap_candi),1,replace=FALSE)
                 swap_one_node<-swap_candi[sample_node]
               }
               prop_incid_matrix<-current_incid_matrix
               prop_incid_matrix[update_order[k],parent_of_update]<-0
               prop_incid_matrix[update_order[k],swap_one_node]<-1
               prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
               prop_trans_func_matrix<-current_trans_func_matrix
               prop_trans_func_matrix[update_order[k],parent_of_update]<-0
               if (aa>0)
                 prop_trans_func_matrix[update_order[k], swap_one_node]<-swap_parent[3]
               if (aa==0)
                 {
                  func_order<-10+sample.int(2,1)
                  prop_trans_func_matrix[update_order[k], swap_one_node]<-func_order
                 }
               swap_one_prob<-1/length(swap_candi)
               prop_trans_func_prob<-1/2  #  new proposal transition function matrix
               xxx<-Error_LLH(prop_trans_func_matrix)
               prop_sample<-xxx[[2]]
               prop_post<-xxx[[1]][length(xxx[[1]])]
               prop_sample_prob<-numeric(); p2c_prob<-numeric()  # calculate acceptance probability
               curr_sample_prob<-numeric(); c2p_prob<-numeric()
               prop_sample_prob[1]<-prop_post;prop_sample_prob[2]<-0       #  this is the prior & likelihood
               curr_sample_prob[1]<-current_post[1];curr_sample_prob[2]<-0
               prop_sample_prob[3]<-log(prop_trans_func_prob); curr_sample_prob[3]<-log(1)    #  this is P(R|T)

               if (length(legal_parent)>0)
                 p2c_prob[1]<-1/3 #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
               if (length(legal_parent)==0)
                 p2c_prob[1]<-1/2

               if (length(legal_parent)>0)
                 c2p_prob[1]<-1/3*swap_one_prob   # this is Q(T_p|T_c) probability of which move and which node selected
               if (length(legal_parent)==0)
                 c2p_prob[1]<-1/2*swap_one_prob
               p2c_prob[2]<-1   # Q(R_c|T_c)
               c2p_prob[2]<-prop_trans_func_prob  # Q(R_p|T_p)
               p2c_prob[3]<-1# Q(\theta_c|T_c, R_c) no root node is introduced.
               c2p_prob[3]<-1                     # Q(\theta_p|T_p, R_p)

               nume<-sum((prop_sample_prob))+sum(log(p2c_prob))
               deno<-sum((curr_sample_prob))+sum(log(c2p_prob))
               acce_prob<-exp(nume-deno); ratio<-runif(1)
               if (ratio<=acce_prob)
                  n<-n+1
             }
           ##################################
           # proposal move 3; remove one parent ### this move may introduce new root node
           if(uu>=2/6 & uu<3/6)
            {
             remove_one_node<-parent_of_update; remove_one_prob<-1
             prop_incid_matrix<-current_incid_matrix
             prop_incid_matrix[update_order[k],remove_one_node]<-0  #remove parent move does not need function.
             prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
             prop_trans_func_matrix<-current_trans_func_matrix
             prop_trans_func_matrix[update_order[k], remove_one_node]<-0
             prop_trans_func_prob<-1
             xxx<-Error_LLH(prop_trans_func_matrix)
             prop_sample<-xxx[[2]]
             prop_post<-xxx[[1]][length(xxx[[1]])]
             prop_sample_prob<-numeric(); p2c_prob<-numeric()   # calculate acceptance probability
             curr_sample_prob<-numeric(); c2p_prob<-numeric()
             prop_sample_prob[1]<-prop_post;prop_sample_prob[2]<-0       #  this is the prior & likelihood
             curr_sample_prob[1]<-current_post[1];curr_sample_prob[2]<-0
             prop_sample_prob[3]<-log(prop_trans_func_prob); curr_sample_prob[3]<-log(1)    #  this is P(R|T)

             if (length(legal_parent)==0)
              p2c_prob[1]<-1
             if (length(legal_parent)==1)
              p2c_prob[1]<-1 #  this is Q(T_c|T_p) inverse move
             if (length(legal_parent)>1)
              p2c_prob[1]<-1/2

             if (length(legal_parent)>0)
               c2p_prob[1]<-1/3*remove_one_prob   # this is Q(T_p|T_c) probability of which move and which node selected
             if (length(legal_parent)==0)
               c2p_prob[1]<-1/2*remove_one_prob
             p2c_prob[2]<-1/2  # Q(R_c|T_c)
             c2p_prob[2]<-prop_trans_func_prob  # Q(R_p|T_p)
             p2c_prob[3]<-1# Q(\theta_c|T_c, R_c)
             c2p_prob[3]<-1       # Q(\theta_p|T_p, R_p)

             nume<-sum((prop_sample_prob))+sum(log(p2c_prob))
             deno<-sum((curr_sample_prob))+sum(log(c2p_prob))
             acce_prob<-exp(nume-deno)
             ratio<-runif(1)
             if (ratio<=acce_prob)
                n<-n+1
           }
           ####################################
           # additional move 1: reverse one arc in pairwise genes
          if(uu>=3/6 & uu<4/6)
            {
             prop_incid_matrix<-current_incid_matrix
             prop_incid_matrix[update_order[k],parent_of_update]<-0
             prop_incid_matrix[parent_of_update,update_order[k]]<-1
             parent_parent<-numeric(); j<-0
             for (i in 1:nrow(current_incid_matrix))
               if (current_incid_matrix[parent_of_update,i]>0)
                 {
                  j<-j+1; parent_parent[j]<-i
                 }
             prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
             if(check.ances.matrix(prop_ances_matrix)==0 & length(parent_parent)==0)  # make sure no loops & parent_of_update no parents
              {
               prop_trans_func_matrix<-current_trans_func_matrix
               func_order<-current_trans_func_matrix[update_order[k],parent_of_update]
               prop_trans_func_matrix[update_order[k],parent_of_update]<-0
               prop_trans_func_matrix[parent_of_update,update_order[k]]<-func_order
               xxx<-Error_LLH(prop_trans_func_matrix)
               prop_sample<-xxx[[2]]
               prop_post<-xxx[[1]][length(xxx[[1]])]

               nume<-sum((prop_post))
               deno<-sum((current_post))
               acce_prob<-exp(nume-deno); ratio<-runif(1)
               if (ratio<=acce_prob)
                  n<-n+1
              }
            }
            ###################################
            # additional move 2: reverse two arcs simultanously
           if(uu>=4/6 & uu<5/6)
            {
             parent_parent<-numeric(); j1<-0
             children<-numeric(); j<-0           # find the children of parent_of_update
             for (i  in 1:ncol(current_incid_matrix))
              {
               if (i !=update_order[k] & current_incid_matrix[i, parent_of_update]>0)
                {
                 j<-j+1; children[j]<-i
                }
               if (current_incid_matrix[parent_of_update, i]>0)
                {
                 j1<-j1+1; parent_parent[j1]<-i
                }
              }
             if (length(children)==1 & length(parent_parent)==0)   # parent_of_update has two exact 2 children and no parents
              {
               parent_child<-numeric(); j<-0        # find parent of children to determine the trans_func_matrix
               for (i in 1:num.node)
                 if (current_incid_matrix[children, i]>0)
                 {
                  j<-j+1; parent_child[j]<-i
                 }

               prop_incid_matrix<-current_incid_matrix
               prop_incid_matrix[update_order[k], parent_of_update]<-0
               prop_incid_matrix[children, parent_of_update]<-0
               prop_incid_matrix[parent_of_update, update_order[k]]<-1
               prop_incid_matrix[parent_of_update, children]<-1
               prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
               if (check.ances.matrix(prop_ances_matrix)==0)  # ensure no directed cycles
               {
                prop_trans_func_matrix<-current_trans_func_matrix
                prop_trans_func_matrix[update_order[k], parent_of_update]<-0
                prop_trans_func_matrix[children, parent_of_update]<-0
                func_order<-sample.int(10,1)
                prop_trans_func_matrix[parent_of_update, update_order[k]]<-func_order
                prop_trans_func_matrix[parent_of_update, children]<-func_order
                if (length(parent_child)==2)
                    prop_trans_func_matrix[children,setdiff(parent_child, parent_of_update)]<-10+sample.int(2,1)
                xxx<-Error_LLH(prop_trans_func_matrix)
                prop_sample<-xxx[[2]]
                prop_post<-xxx[[1]][length(xxx[[1]])]

                nume<-sum((prop_post))
                deno<-sum((current_post))
                acce_prob<-exp(nume-deno); ratio<-runif(1)
                if (ratio<=acce_prob)
                  n<-n+1
               }
             }
           }
         ############################################# additional move: swap two with one   
           if (uu>=5/6 & length(swap_candi)>1)
             {
              legal.triplet.pare<-intersect(swap_candi, triplet.prior.pare)
              prop.prob<-runif(1); aa<-0
              if(prop.prob>prop.ratio)
               {
                swap.candidate<-list()
                if (length(intersect(parent_of_update, legal.triplet.pare))==0 & nrow(triplet.prior.set)>0)
                 for (jj in 1: nrow(triplet.prior.pare))
                  if (length(intersect(parent_of_update, triplet.prior.pare[jj,]))==0 & length(intersect(triplet.prior.pare[jj,], swap_candi))==2)
                   {
                    aa<-aa+1
                    swap.candidate[[aa]]<-triplet.prior.set[jj,]
                   }
                 if (length(swap.candidate)>0)
                   {
                     score<-numeric()
                     for (jj in 1:length(swap.candidate))
                       score[jj]<-swap.candidate[[jj]][5]
                     bbb<-matrix(nrow=length(swap.candidate),ncol=6) # convert list to matrix
                     for (jjj in 1:length(swap.candidate))
                       bbb[jjj,]<-swap.candidate[[jjj]]
                     max.score.candidate<-bbb[bbb[,5]==max(score),]
                     if (is.numeric(max.score.candidate)==T)
                       swap_parent<-max.score.candidate
                     if (is.matrix(max.score.candidate)==T)
                      if (nrow(max.score.candidate)>1)
                       swap_parent<-max.score.candidate[sample.int(nrow(max.score.candidate),1),]
                     swap_two_parent<-swap_parent[1:2]
                   }
                }
              if (prop.prob<prop.ratio||aa==0)
                {
                  sample_two_node<-sample.int(length(swap_candi),2,replace=FALSE)
                  swap_two_parent<-c(swap_candi[sample_two_node[1]], swap_candi[sample_two_node[2]])
                }
              prop_incid_matrix<-current_incid_matrix
              prop_incid_matrix[update_order[k],parent_of_update]<-0
              prop_incid_matrix[update_order[k],swap_two_parent[1]]<-1
              prop_incid_matrix[update_order[k],swap_two_parent[2]]<-1
              prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
              prop_trans_func_matrix<-current_trans_func_matrix
              prop_trans_func_matrix[update_order[k],parent_of_update]<-0
              if (aa>0)
                {
                 prop_trans_func_matrix[update_order[k], swap_two_parent[1]]<-swap_parent[4]
                 prop_trans_func_matrix[update_order[k], swap_two_parent[2]]<-swap_parent[4]
                }
              if (aa==0)
                {
                 func_order<-sample.int(10,1)
                 prop_trans_func_matrix[update_order[k], swap_two_parent[1]]<-func_order
                 prop_trans_func_matrix[update_order[k], swap_two_parent[2]]<-func_order
                }
              swap_two_prob<-1/choose(length(swap_candi),2)
              prop_trans_func_prob<-1/10
              xxx<-Error_LLH(prop_trans_func_matrix)
              prop_sample<-xxx[[2]]
              prop_post<-xxx[[1]][length(xxx[[1]])]
              prop_sample_prob<-numeric(); p2c_prob<-numeric() # calculate acceptance probability
              curr_sample_prob<-numeric(); c2p_prob<-numeric()

              prop_sample_prob[1]<-prop_post;prop_sample_prob[2]<-0       #  this is the prior & likelihood
              curr_sample_prob[1]<-current_post[1];curr_sample_prob[2]<-0
              prop_sample_prob[3]<-log(prop_trans_func_prob); curr_sample_prob[3]<-log(1)    #  this is P(R|T)

              if (length(swap_candi)>1)
                p2c_prob[1]<-1/4 #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
              if (length(swap_candi)==1)
                p2c_prob[1]<-1/3

              if (length(swap_candi)>1)
                c2p_prob[1]<-1/4*swap_two_prob   # this is Q(T_p|T_c) probability of which move and which node selected
              if (length(swap_candi)==1)
                c2p_prob[1]<-1/3*swap_two_prob
              p2c_prob[2]<-1/10   # Q(R_c|T_c)
              c2p_prob[2]<-prop_trans_func_prob  # Q(R_p|T_p)
              p2c_prob[3]<-1# Q(\theta_c|T_c, R_c) no root node is introduced.
              c2p_prob[3]<-1                     # Q(\theta_p|T_p, R_p)

              nume<-sum((prop_sample_prob))+sum(log(p2c_prob))
              deno<-sum((curr_sample_prob))+sum(log(c2p_prob))
              acce_prob<-exp(nume-deno); ratio<-runif(1)
              if (ratio<=acce_prob)
                 n<-n+1
             }
        }
#########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# case3: two parents
       if (length(parent_of_update)==2)
         {
           uu<-runif(1); #two.parent.ratio=2/7*(1/BIC.pairwise/(1/BIC.pairwise+1/BIC.standalone))
           two.parent.ratio=1/7
           ##########################
           #proposal move 1: remove one parent
          if (uu<two.parent.ratio)
            {
                 pare.prior.overlap<-intersect(parent_of_update, pairwise.prior.pare)
                 prop.prob<-runif(1);aa<-0
                 if(prop.prob>=prop.ratio)
                  if (length(pare.prior.overlap)==1 & nrow(pairwise.prior.set)>0)
                  {
                   aa<-aa+1
                   remove_one_node<-setdiff(parent_of_update, pairwise.prior.pare)
                   remove_pare<-pairwise.prior.set[pairwise.prior.set[,1]==setdiff(parent_of_update,remove_one_node),]
                   if (is.numeric(remove_pare)==T)
                     func_order<-remove_pare[3]
                   if(is.matrix(remove_pare)==T)
                    {
                     Remove_Pare<-remove_pare[remove_pare[,4]==max(remove_pare[,4]),]
                     if (nrow(Remove_Pare)>1)
                       Remove_Pare<-Remove_Pare[sample.int(nrow(Remove_Pare), 1),]  # each candidate has equal chance
                     func_order<-Remove_Pare[3]
                    }
                   }
                 if (prop.prob<prop.ratio||aa==0)
                   {
                    sample_node<-sample.int(length(parent_of_update),1,replace=FALSE)
                    remove_one_node<-parent_of_update[sample_node]
                   }
                 remove_one_prob<-1/length(parent_of_update)
                 prop_incid_matrix<-current_incid_matrix; prop_incid_matrix[update_order[k], remove_one_node]<-0
                 prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
                 prop_trans_func_matrix<-current_trans_func_matrix
                 prop_trans_func_matrix[update_order[k], remove_one_node]<-0
                 if (aa>0)
                   prop_trans_func_matrix[update_order[k], setdiff(parent_of_update, remove_one_node)]<-func_order # by the definition
                 if (aa==0)
                    func_order<-10+sample.int(2,1)
                 prop_trans_func_matrix[update_order[k], setdiff(parent_of_update, remove_one_node)]<-func_order
                 prop_trans_func_prob<-1/2
                 xxx<-Error_LLH(prop_trans_func_matrix)
                 prop_sample<-xxx[[2]]
                 prop_post<-xxx[[1]][length(xxx[[1]])]
                 prop_sample_prob<-numeric(); p2c_prob<-numeric() # calculate acceptance probability
                 curr_sample_prob<-numeric(); c2p_prob<-numeric()

                 prop_sample_prob[1]<-prop_post;prop_sample_prob[2]<-0       #  this is the prior & likelihood
                 curr_sample_prob[1]<-current_post[1];curr_sample_prob[2]<-0
                 prop_sample_prob[3]<-log(prop_trans_func_prob); curr_sample_prob[3]<-log(1)    #  this is P(R|T)

                 if (length(swap_candi)>0)
                  p2c_prob[1]<-1/3 #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
                 if (length(swap_candi)==0)
                  p2c_prob[1]<-1/2

                 if (length(swap_candi)==0)
                  c2p_prob[1]<-1/2
                 if (length(swap_candi)>1)
                  c2p_prob[1]<-1/4*remove_one_prob   # this is Q(T_p|T_c) probability of which move and which node selected
                 if (length(swap_candi)==1)
                   c2p_prob[1]<-1/3*remove_one_prob

                 p2c_prob[2]<-1/10   # Q(R_c|T_c)
                 c2p_prob[2]<-prop_trans_func_prob  # Q(R_p|T_p)
                 p2c_prob[3]<-1  # Q(\theta_c|T_c, R_c) no root node is introduced.
                 c2p_prob[3]<-1                     # Q(\theta_p|T_p, R_p)

                 nume<-sum((prop_sample_prob))+sum(log(p2c_prob))
                 deno<-sum((curr_sample_prob))+sum(log(c2p_prob))
                 acce_prob<-exp(nume-deno)
                 ratio<-runif(1)
                 if (ratio<=acce_prob)
                    n<-n+1
             }
          ############################
          # proposal move 2; remove two parents
           if (uu>=two.parent.ratio & uu<2/7)
            {
             remove_node<-parent_of_update;  remove_prob<-1
             prop_incid_matrix<-current_incid_matrix
             prop_incid_matrix[update_order[k], remove_node[1]]<-0; prop_incid_matrix[update_order[k], remove_node[2]]<-0
             prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
             prop_trans_func_matrix<-current_trans_func_matrix
             prop_trans_func_matrix[update_order[k], remove_node[1]]<-0; prop_trans_func_matrix[update_order[k], remove_node[2]]<-0
             prop_trans_func_prob<-1
             xxx<-Error_LLH(prop_trans_func_matrix)
             prop_sample<-xxx[[2]]
             prop_post<-xxx[[1]][length(xxx[[1]])]
             prop_sample_prob<-numeric(); p2c_prob<-numeric()    # calculate acceptance probability
             curr_sample_prob<-numeric(); c2p_prob<-numeric()
             prop_sample_prob[1]<-prop_post;prop_sample_prob[2]<-0       #  this is the prior & likelihood
             curr_sample_prob[1]<-current_post[1];curr_sample_prob[2]<-0
             prop_sample_prob[3]<-log(prop_trans_func_prob); curr_sample_prob[3]<-log(1)    #  this is P(R|T)

             p2c_prob[1]<-1/2 #  this is Q(T_c|T_p) inverse move
             c2p_prob[1]<-1/4*remove_prob   # this is Q(T_p|T_c) probability of which move and which node selected
             p2c_prob[2]<-1/10   # Q(R_c|T_c)
             c2p_prob[2]<-prop_trans_func_prob  # Q(R_p|T_p)
             p2c_prob[3]<-1# Q(\theta_c|T_c, R_c) no root node is introduced.
             c2p_prob[3]<-1       # Q(\theta_p|T_p, R_p)

             nume<-sum((prop_sample_prob))+sum(log(p2c_prob))
             deno<-sum((curr_sample_prob))+sum(log(c2p_prob))
             acce_prob<-exp(nume-deno); ratio<-runif(1)
             if (ratio<=acce_prob)
                n<-n+1
           }
          ############################
          #proposal move 3: swap one parent
            if (uu>=2/7 & uu<3/7 & length(swap_candi)>0)
            {
              legal.triplet.pare<-intersect(swap_candi, triplet.prior.pare)
              prop.prob<-runif(1);aa<-0
              if(prop.prob>=prop.ratio)
               {
                swap.candidate<-list()
                if (length(intersect(parent_of_update, legal.triplet.pare))==1 & nrow(triplet.prior.set)>0)
                 for (jj in 1: nrow(triplet.prior.pare))
                  if (length(intersect(parent_of_update, triplet.prior.pare[jj,]))==1)
                   {
                    aa<-aa+1
                    swap.candidate[[aa]]<-triplet.prior.set[jj,]
                    }
                 if (length(swap.candidate)>0)
                 {
                  score<-numeric()
                  for (jj in 1:length(swap.candidate))
                    score[jj]<-swap.candidate[[jj]][5]
                  bbb<-matrix(nrow=length(swap.candidate),ncol=6)
                  for (jjj in 1:length(swap.candidate))
                    bbb[jjj,]<-swap.candidate[[jjj]]

                  max.score.candidate<-bbb[bbb[,5]==max(score),]
                  if (is.numeric(max.score.candidate)==T)
                    swap_parent<-max.score.candidate
                  if (is.matrix(max.score.candiate)==T)
                   if (nrow(max.score.candidate)>1)
                    swap_parent<-max.score.candidate[sample.int(nrow(max.score.candidate),1),]
                  swap_one_node<-setdiff(swap_parent[1:2], parent_of_update)
                  sample_one_parent<-setdiff(parent_of_update,swap_parent[1:2])
                 }
                }
              if (prop.prob<prop.ratio||aa==0)
               {
                sample_node<-sample.int(length(swap_candi),1,replace=FALSE); sample_parent<-sample.int(2,1,replace=FALSE)
                swap_one_node<-swap_candi[sample_node]; sample_one_parent<-parent_of_update[sample_parent]
               }
               prop_incid_matrix<-current_incid_matrix
               prop_incid_matrix[update_order[k],sample_one_parent]<-0
               prop_incid_matrix[update_order[k],swap_one_node]<-1
               prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
               prop_trans_func_matrix<-current_trans_func_matrix
               prop_trans_func_matrix[update_order[k],sample_one_parent]<-0
               if (aa>0)
                {
                 prop_trans_func_matrix[update_order[k], setdiff(parent_of_update, sample_one_parent)]<-swap_parent[4]
                 prop_trans_func_matrix[update_order[k], swap_one_node]<-swap_parent[4]
                }
               if (aa==0)
                {
                 func_order<-sample.int(10,1)
                 prop_trans_func_matrix[update_order[k], setdiff(parent_of_update, sample_one_parent)]<-func_order
                 prop_trans_func_matrix[update_order[k], swap_one_node]<-func_order
                }
               swap_one_prob<-1/length(swap_candi)*1/2
               prop_trans_func_prob<-1/10
               xxx<-Error_LLH(prop_trans_func_matrix)
               prop_sample<-xxx[[2]]
               prop_post<-xxx[[1]][length(xxx[[1]])]
               prop_sample_prob<-numeric(); p2c_prob<-numeric() # calculate acceptance probability
               curr_sample_prob<-numeric(); c2p_prob<-numeric()
               prop_sample_prob[1]<-prop_post;prop_sample_prob[2]<-0       #  this is the prior & likelihood
               curr_sample_prob[1]<-current_post[1];curr_sample_prob[2]<-0
               prop_sample_prob[3]<-log(prop_trans_func_prob); curr_sample_prob[3]<-log(1)    #  this is P(R|T)

               if (length(swap_candi)==0)
                 p2c_prob[1]<-1   # since we use log(p2c_prob)
               if (length(swap_candi)>1)
                 p2c_prob[1]<-1/4 #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
               if (length(swap_candi)==1)
                 p2c_prob[1]<-1/3

               if (length(swap_candi)==0)
                 c2p_prob[1]<-1
               if (length(swap_candi)>1)
                 c2p_prob[1]<-1/4*swap_one_prob   # this is Q(T_p|T_c) probability of which move and which node selected
               if (length(swap_candi)==1)
                 c2p_prob[1]<-1/3*swap_one_prob
               p2c_prob[2]<-1/10   # Q(R_c|T_c)
               c2p_prob[2]<-prop_trans_func_prob  # Q(R_p|T_p)
               p2c_prob[3]<-1# Q(\theta_c|T_c, R_c) no root node is introduced.
               c2p_prob[3]<-1                     # Q(\theta_p|T_p, R_p)

               nume<-sum((prop_sample_prob))+sum(log(p2c_prob))
               deno<-sum((curr_sample_prob))+sum(log(c2p_prob))
               acce_prob<-exp(nume-deno); ratio<-runif(1)
               if (ratio<=acce_prob)
                 n<-n+1
            }
         ##############################
         #proposal move 4, swap two parents
           if (uu>=3/7 & uu<4/7 & length(swap_candi)>1)
           {
              legal.triplet.pare<-intersect(swap_candi, triplet.prior.pare)
              prop.prob<-runif(1); aa<-0
              if(prop.prob>prop.ratio)
               {
                swap.candidate<-list()
                if (length(intersect(parent_of_update, legal.triplet.pare))==0 & nrow(triplet.prior.set)>0)
                 for (jj in 1: nrow(triplet.prior.pare))
                  if (length(intersect(parent_of_update, triplet.prior.pare[jj,]))==0 & length(intersect(triplet.prior.pare[jj,], swap_candi))==2)
                   {
                    aa<-aa+1
                    swap.candidate[[aa]]<-triplet.prior.set[jj,]
                   }
                 if ( length(swap.candidate)>0)
                   {
                     score<-numeric()
                     for (jj in 1:length(swap.candidate))
                       score[jj]<-swap.candidate[[jj]][5]
                     bbb<-matrix(nrow=length(swap.candidate),ncol=6) # convert list to matrix
                     for (jjj in 1:length(swap.candidate))
                       bbb[jjj,]<-swap.candidate[[jjj]]
                     max.score.candidate<-bbb[bbb[,5]==max(score),]
                     if (is.numeric(max.score.candidate)==T)
                       swap_parent<-max.score.candidate
                     if (is.matrix(max.score.candidate)==T)
                      if (nrow(max.score.candidate)>1)
                       swap_parent<-max.score.candidate[sample.int(nrow(max.score.candidate),1),]
                     swap_two_parent<-swap_parent[1:2]
                   }
                }
              if (prop.prob<prop.ratio||aa==0)
                {
                  sample_two_node<-sample.int(length(swap_candi),2,replace=FALSE)
                  swap_two_parent<-c(swap_candi[sample_two_node[1]], swap_candi[sample_two_node[2]])
                }
              prop_incid_matrix<-current_incid_matrix
              prop_incid_matrix[update_order[k],parent_of_update[1]]<-0
              prop_incid_matrix[update_order[k],parent_of_update[2]]<-0
              prop_incid_matrix[update_order[k],swap_two_parent[1]]<-1
              prop_incid_matrix[update_order[k],swap_two_parent[2]]<-1
              prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
              prop_trans_func_matrix<-current_trans_func_matrix
              prop_trans_func_matrix[update_order[k],parent_of_update[1]]<-0
              prop_trans_func_matrix[update_order[k],parent_of_update[2]]<-0
              if (aa>0)
                {
                 prop_trans_func_matrix[update_order[k], swap_two_parent[1]]<-swap_parent[4]
                 prop_trans_func_matrix[update_order[k], swap_two_parent[2]]<-swap_parent[4]
                }
              if (aa==0)
                {
                 func_order<-sample.int(10,1)
                 prop_trans_func_matrix[update_order[k], swap_two_parent[1]]<-func_order
                 prop_trans_func_matrix[update_order[k], swap_two_parent[2]]<-func_order
                }
              swap_two_prob<-1/choose(length(swap_candi),2)
              prop_trans_func_prob<-1/10
              xxx<-Error_LLH(prop_trans_func_matrix)
              prop_sample<-xxx[[2]]
              prop_post<-xxx[[1]][length(xxx[[1]])]
              prop_sample_prob<-numeric(); p2c_prob<-numeric() # calculate acceptance probability
              curr_sample_prob<-numeric(); c2p_prob<-numeric()

              prop_sample_prob[1]<-prop_post;prop_sample_prob[2]<-0       #  this is the prior & likelihood
              curr_sample_prob[1]<-current_post[1];curr_sample_prob[2]<-0
              prop_sample_prob[3]<-log(prop_trans_func_prob); curr_sample_prob[3]<-log(1)    #  this is P(R|T)

              if (length(swap_candi)>1)
                p2c_prob[1]<-1/4 #  this is Q(T_c|T_p) inverse move, because there are 4 moves after adding one parents, depending on the number of nodes to swap
              if (length(swap_candi)==1)
                p2c_prob[1]<-1/3

              if (length(swap_candi)>1)
                c2p_prob[1]<-1/4*swap_two_prob   # this is Q(T_p|T_c) probability of which move and which node selected
              if (length(swap_candi)==1)
                c2p_prob[1]<-1/3*swap_two_prob
              p2c_prob[2]<-1/10   # Q(R_c|T_c)
              c2p_prob[2]<-prop_trans_func_prob  # Q(R_p|T_p)
              p2c_prob[3]<-1# Q(\theta_c|T_c, R_c) no root node is introduced.
              c2p_prob[3]<-1                     # Q(\theta_p|T_p, R_p)

              nume<-sum((prop_sample_prob))+sum(log(p2c_prob))
              deno<-sum((curr_sample_prob))+sum(log(c2p_prob))
              acce_prob<-exp(nume-deno); ratio<-runif(1)
              if (ratio<=acce_prob)
                 n<-n+1
            }
        #############################################################
       #  additional move 1: reorder the input and output variables  for XOR relation
             if (uu>=4/7 & uu<5/7 & current_trans_func_matrix[update_order[k],parent_of_update[1]]==9)
             {
              prop_incid_matrix<-current_incid_matrix
              sample.child<-parent_of_update[sample.int(2,1)]
              for (i in 1:2)
                prop_incid_matrix[update_order[k],parent_of_update[i]]<-0
              parent_parent<-numeric(); mm<-0
              for (i in 1:num.node)
                if (current_incid_matrix[sample.child,i ]>0)
                   {
                    mm<-mm+1; parent_parent[mm]<-i
                   }
             if (length(parent_parent)==0)      # selected child should have no parents before this operation
              {
              prop_incid_matrix[sample.child,update_order[k]]<-1
              prop_incid_matrix[sample.child,setdiff(parent_of_update,sample.child)]<-1
              prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
              if (check.ances.matrix(prop_ances_matrix)==0)      # make sure no directed cycles
               {
                prop_trans_func_matrix<-current_trans_func_matrix
                func_order<-current_trans_func_matrix[update_order[k],parent_of_update[1]]
                for (i in 1:2)
                  prop_trans_func_matrix[update_order[k],parent_of_update[i]]<-0
                prop_trans_func_matrix[sample.child,update_order[k]]<-func_order
                prop_trans_func_matrix[sample.child,setdiff(parent_of_update,sample.child)]<-func_order
                xxx<-Error_LLH(prop_trans_func_matrix)
                prop_sample<-xxx[[2]]
                prop_post<-xxx[[1]][length(xxx[[1]])]

                nume<-sum((prop_post))
                deno<-sum((current_post))
                acce_prob<-exp(nume-deno); ratio<-runif(1)
                if (ratio<=acce_prob)
                   n<-n+1
                }
               }
             }
       ###################################################
       # additional move 2: reverse one arc among gene triplet
         if (uu>=5/7 & uu<6/7)
             {
              reverse_pare<-parent_of_update[sample.int(2,1)]
              remain_pare<-setdiff(parent_of_update, reverse_pare)
               ###################
              reverse_node_pare<-numeric();j<-1         # find the already existing parent for reverse_pare
              for (i in 1: ncol(current_incid_matrix))
                if (current_incid_matrix[reverse_pare,i]!=0)
                 {
                  reverse_node_pare[j]<-i; j<-j+1
                 }
              if (length(reverse_node_pare)==0)    # currently no parent
               {
                prop_incid_matrix<-current_incid_matrix
                prop_incid_matrix[update_order[k], reverse_pare]<-0
                prop_incid_matrix[reverse_pare, update_order[k]]<-1
                prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
                if (check.ances.matrix(prop_ances_matrix)==0) # make sure no loops
                 {
                  prop_trans_func_matrix<-current_trans_func_matrix
                  prop_trans_func_matrix[update_order[k], reverse_pare]<-0
                  ###############################
                  prop_trans_func_matrix[update_order[k], remain_pare]<-10+sample.int(2,1)   # try no use prior
                  prop_trans_func_matrix[reverse_pare, update_order[k]]<-10+sample.int(2,1)
                  xxx<-Error_LLH(prop_trans_func_matrix)
                  prop_sample<-xxx[[2]]
                  prop_post<-xxx[[1]][length(xxx[[1]])]

                  nume<-sum((prop_post))
                  deno<-sum((current_post))
                  acce_prob<-exp(nume-deno); ratio<-runif(1)
                  if (ratio<=acce_prob)
                    n<-n+1
                 }
              }
             if (length(reverse_node_pare)==1)       # currently one parent
              {
               prop_incid_matrix<-current_incid_matrix
               prop_incid_matrix[update_order[k], reverse_pare]<-0
               prop_incid_matrix[reverse_pare, update_order[k]]<-1
               prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
               if (check.ances.matrix(prop_ances_matrix)==0) # make sure no loops
               {
                prop_trans_func_matrix<-current_trans_func_matrix
                prop_trans_func_matrix[update_order[k], reverse_pare]<-0
                ###############################
                prop_trans_func_matrix[update_order[k], remain_pare]<-10+sample.int(2,1)   # no use prior
                func_order<-sample.int(10,1)
                prop_trans_func_matrix[reverse_pare, update_order[k]]<-func_order
                prop_trans_func_matrix[reverse_pare, reverse_node_pare]<-func_order
                xxx<-Error_LLH(prop_trans_func_matrix)
                prop_sample<-xxx[[2]]
                prop_post<-xxx[[1]][length(xxx[[1]])]

                nume<-sum((prop_post))
                deno<-sum((current_post))
                acce_prob<-exp(nume-deno); ratio<-runif(1)
                if (ratio<=acce_prob)
                  n<-n+1
               }
              }
          }
       ########################################
       # additional move 3: only change boolean functions among gene triplets  note: the network may remain the same
            if (uu>=6/7)
             {
              prop.prob<-runif(1);aa<-0
              if(prop.prob>=prop.ratio)
                if (nrow(triplet.prior.set)>=1 & ncol(triplet.prior.set)>1)
                  {
                   candidate.pare.set<-list()
                   for (jj in 1:nrow(triplet.prior.set))
                     if (length(intersect(parent_of_update,triplet.prior.pare[jj,] ))==2)
                      {
                       aa<-aa+1
                       candidate.pare.set[[aa]]<-triplet.prior.set[jj,]
                      }
                   if (length(candidate.pare.set)>0)
                    {
                     score<-numeric()
                     for (jj in 1:length(candidate.pare.set))
                      score[jj]<-candidate.pare.set[[jj]][5]
                     bbb<-matrix(nrow=length(candidate.pare.set),ncol=6)
                     for (jjj in 1:length(candidate.pare.set))
                      bbb[jjj,]<-candidate.pare.set[[jjj]]

                     max.score.candidate<-bbb[bbb[,5]==max(score),]
                     if (is.numeric(max.score.candidate)==T)
                       add_parent<-max.score.candidate
                     if (is.matrix(max.score.candidate)==T)
                      if (nrow(max.score.candidate)>1)
                       add_parent<-max.score.candidate[sample.int(nrow(max.score.candidate),1),]
                    }
                  }
              prop_incid_matrix<-current_incid_matrix
              prop_ances_matrix<-update.ancestor_matrix(prop_incid_matrix)
              prop_trans_func_matrix<-current_trans_func_matrix
              if (aa>0)
                for (i in 1:2)
                  prop_trans_func_matrix[update_order[k], parent_of_update[i]]<-add_parent[4]
              if (prop.prob<prop.ratio||aa==0)
               {
                func_order<-sample.int(10,1)
                for (i in 1:2)
                  prop_trans_func_matrix[update_order[k], parent_of_update[i]]<-func_order
               }
              xxx<-Error_LLH(prop_trans_func_matrix)
              prop_sample<-xxx[[2]]
              prop_post<-xxx[[1]][length(xxx[[1]])]
              nume<-sum(prop_post)
              deno<-sum(current_post)
              acce_prob<-exp(nume-deno); ratio<-runif(1)
              if ( acce_prob!=1 & ratio<=acce_prob)
                n<-n+1
             }
        }
#######################################
#######################################
       if (old==n)    # no update
        {
         iter<-iter+1
         xxx<-Error_LLH(current_trans_func_matrix)
         current_sample<-xxx[[2]]
         Sample_Matrix[iter,]<-current_sample
         all_logpost[iter]<-xxx[[1]][length(xxx[[1]])]
         Incidence_Matrix[[iter]]<-current_incid_matrix
         Ancestor_Matrix[[iter]]<-current_ances_matrix
         Trans_Func_Matrix[[iter]]<-current_trans_func_matrix
        }
       if (old!=n)  # when either T or F changes,  recount the number
        {
         iter<-iter+1
         jump_point[n]<-iter

         logpost[n]<-sum(prop_post)
         Sample_Matrix[iter,]<-prop_sample
         all_logpost[iter]<-sum(prop_post)
         Incidence_Matrix[[iter]]<-prop_incid_matrix
         Ancestor_Matrix[[iter]]<-prop_ances_matrix
         Trans_Func_Matrix[[iter]]<-prop_trans_func_matrix
        }
    }  # end of updating

 } # end of num_update
##########################
#num<-numeric()  # summarize the result
#for (i in 2:length(jump_point))
#     num[i-1]<-jump_point[i]-jump_point[i-1]
#num[length(jump_point)]<-iter+1-jump_point[length(jump_point)]
#for (i in 1:length(num))
#  if (num[i]==max(num))
#     max_num<-i
######################################
true_root<-numeric(); j<-0
for ( i in 1:nrow(true_network))
  if (sum(true_network[i,])==0)
   {
    j<-j+1; true_root[j]<-i
   }
common_links<-list(); common_root<-list(); model_root<-list()
non_inferred_links<-numeric(); distance_matrix1<-list(); false_links<-numeric()
correct_links_ratio1<-numeric(); correct_links_ratio2<-numeric()
for (i in 1:length(Incidence_Matrix))
  {
   k<-0 ; root_node<-numeric()
   for (j in 1:nrow(Incidence_Matrix[[i]]))
     if (sum(Incidence_Matrix[[i]][j,])==0)
       {
        k<-k+1; root_node[k]<-j
       }
   model_root[[i]]<-root_node
   common_root[[i]]<-intersect(true_root, root_node)
   distance_matrix1[[i]]<-true_incid_matrix-Incidence_Matrix[[i]]
   bb<-distance_matrix1[[i]]
   cc<-bb[bb>0]; dd<-bb[bb<0]
   non_inferred_links[i]<-sum(cc)
   false_links[i]<-abs(sum(dd))
   correct_links_ratio1[i]<-(sum(Incidence_Matrix[[i]])-false_links[i]+length(common_root[[i]]))/(sum(Incidence_Matrix[[i]])+length(root_node))  # correct_links_ratio=correctly inferred links/inferred total links
   correct_links_ratio2[i]<-(sum(true_incid_matrix)-non_inferred_links[i]+length(common_root[[i]]))/(sum(true_incid_matrix)+length(true_root))
  }
###########################
total_func<-0        # report inferred structure and functions
for (i in 1:nrow(true_incid_matrix))
  if (sum(true_incid_matrix[i,])>0)
    total_func<-total_func+1
non_inferred_SF<-numeric(); SF_ratio<-numeric();  distance_matrix2<-list()
for (i in 1:length(Trans_Func_Matrix))
  {
   distance_matrix2[[i]]<-true_network-Trans_Func_Matrix[[i]]
   tt<-0
   for (j in 1:nrow(distance_matrix2[[i]]))
     {
      bb<-distance_matrix2[[i]][j,]
      if(max(bb)>0)
        tt<-tt+1
     }
   non_inferred_SF[i]<-tt
   SF_ratio[i]<-(total_func-tt+length(common_root[[i]]))/(total_func+length(true_root))
  }
#  All_Trans_Func_Matrix[[iii]]=Trans_Func_Matrix
#  All_Logpost[[iii]]=all_logpost
#  All_Correct_Rate[[iii]]=SF_ratio
} # end of iii  
time.spent<-proc.time()-ptm
################### test the method
#save(list=ls(),file="BICNoPerfectProposa23.RData")