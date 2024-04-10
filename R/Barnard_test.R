Barnard_test<-function(mat,alternative=c("two.sided","greater","less"),n_p=101L,trace=F,
                       return_data=F,data_previous=NULL){
  t0<-Sys.time()
  mode(mat)<-"integer"
  n1<-sum(mat[1,]);n2<-sum(mat[2,])
  stopifnot(n1>0&n2>0)
  c1<-mat[1,1];c2<-mat[2,1]
  alternative<-match.arg(alternative)
  stopifnot(n_p>=1L)
  if(c1/n1==c2/n2&alternative=="two.sided"){
    if(return_data){
      return(list(p.value=1,data=NULL))  
    }else{
      return(1)
    }
  }
  if(c1/n1>c2/n2){
    if(alternative=="less"){
      c1<-n1-c1
      c2<-n2-c2
    }
  }else{
    if(alternative!="greater"){
      c1<-n1-c1
      c2<-n2-c2
    }
  }
  
  p0<-seq(0,1,length.out=n_p)
  n_p_h<-ceiling(n_p/2)
  n_pm1<-n_p-1L;n11<-n1+1L;n21<-n2+1L
  n0<-n11*n21;n1121<-n11+n21
  
  if(is.null(data_previous)){
    prob_list<-get_prob_list(n_p_h,p0,n1,n11,n2,n21)
  }else{
    stopifnot(data_previous$n1==n1)
    stopifnot(data_previous$n2==n2)
    stopifnot(identical(data_previous$p0,p0))
    stopifnot(data_previous$alternative==alternative)
    prob_list<-data_previous$prob_list
  }
  if(trace){
    cat("Preparing prob_list: done.",capture.output(Sys.time()-t0),"\r\n")
  }
  
  twoside<-alternative=="two.sided"
  select_i0<-c1;select_j0<-c2
  
  if(is.null(data_previous)){
    k<-0L;p_l<-1L
    p_prob<-p_prob_indiv<-vector("list",n1121)
    p_i<-p_j<-rep(NA_integer_,n1121)
    p_prob_max<-rep(NA_real_,n1121)
    p_prob[[1L]]<-p_prob_indiv[[1L]]<-BNDtest_get_prob(n1,0L,n1,n2,n_p,n_pm1, 
                                                       n_p_h,prob_list,twoside)
    p_prob_max[1L]<-max(p_prob[[1L]])
    p_i[1L]<-n1
    p_j[1L]<-0L
    select_i<-integer(0)
    select_j<-integer(0)
  }else{
    select_i<-c(data_previous$select_i)
    select_j<-c(data_previous$select_j)
    loc<-which_loc(select_i0,select_j0,select_i,select_j)
    if(is.na(loc)){
      k<-c(data_previous$k)
      p_l<-c(data_previous$p_l)
      n_each<-length(data_previous$p_prob_indiv[[1]])
      p_prob_indiv<-copy_list(data_previous$p_prob_indiv,p_l,n_each)
      p_prob<-copy_list(data_previous$p_prob,p_l,n_each)
      p_prob_max<-c(data_previous$p_prob_max)
      p_i<-c(data_previous$p_i)
      p_j<-c(data_previous$p_j)
    }else{
      p.value<-max(path_get_prob(select_i,select_j,loc,n1,n2,n_p,n_pm1,n_p_h,prob_list,twoside,trace))
      if(trace){
        cat("\r\n")
        cat("Calculating P value: done.",capture.output(Sys.time()-t0),"\r\n")
      }
      if(return_data){
        return(list(p.value=p.value,data=NULL))  
      }else{
        return(p.value)
      }
    }
  }
  
  out<-BNDtest_loop(p_prob_max,p_l,p_i,p_j,k,n0,trace,select_i0,
                    select_j0,p_prob,p_prob_indiv,n1,n2,select_i,select_j,
                    twoside,n_p,n_pm1,n_p_h,prob_list,p0,alternative)
  
  if(trace){
    cat("\r\n")
    cat("Calculating P value: done.",capture.output(Sys.time()-t0),"\r\n")
  }
  
  if(return_data){
    return(out) 
  }else{
    return(out$p.value)
  }
}
