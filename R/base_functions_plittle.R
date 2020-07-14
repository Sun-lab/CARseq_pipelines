# base functions
smart_mkdir = function(input_dir){
  if( !file.exists(input_dir) || !dir.exists(input_dir) ){
    dir.create(input_dir)
  }
}
smart_RT = function(...){
  #' @title smart_RT
  #' @description Runs \code{\link[utils]{read.table}} but sets \code{stringsAsFactors = FALSE}.
  #' @inherit utils::read.table
  #' @export
  read.table(...,stringsAsFactors=FALSE)
}
smart_WT = function(...){
	#' @title smart_WT
	#' @description Runs \code{\link[utils]{write.table}} but sets \code{row.names = FALSE} and \code{quote = FALSE}.
	#' @inherit utils::write.table
	#' @export
	write.table(...,row.names=FALSE,quote=FALSE)
}
name_change = function(DATA,ORIG_NAME,NEW_NAME){
  new_index = which(names(DATA) == NEW_NAME)
  if( length(new_index) > 0 ){
    DATA = DATA[,-new_index]
  }
  index = which(names(DATA) == ORIG_NAME)
  names(DATA)[index] = NEW_NAME
  DATA
}
smart_merge = function(x,y,mess=NULL,...){
  if( !is.null(mess) ){
    intersect_vars = paste(intersect(names(x),names(y)),collapse=", ")
    cat(paste0("Merging dataframes on variables = { ",intersect_vars," }\n"))
  }
  
  merge(x,y,by=intersect(names(x),names(y)),...)
}
smart_table = function(...){
  table(...,useNA='ifany')
}
smart_df = function(...){
  #' @title smart_df
  #' @description Runs \code{\link[base]{data.frame}} but sets \code{stringsAsFactors = FALSE}.
  #' @inherit base::data.frame
  #' @export
  data.frame(...,stringsAsFactors=FALSE)
}
show_screeplot = function(pca_output,main){
  def_par = par()
  par(mfrow = c(2,1),mar = c(4,4,0,1),oma = c(0,0,2,0))
  barplot(pca_output$values[1:50]/sum(pca_output$values),
          main = "",xlab = "Index",ylab = "Eigen-value")
  barplot(cumsum(pca_output$values)[1:50]/sum(pca_output$values),
          ylim = c(0,1),yaxs = "i",xlab = "Index",
          ylab = "Cumulative Prop. Var. Explained")
  mtext(main,outer = TRUE,cex = 1.5)
  par(mfrow = def_par$mfrow,mar = def_par$mar,oma = def_par$oma)
}
show_pc_color = function(PCS,DATA = NULL,VAR = NULL,CAT,submain){
  if( is.null(VAR) ){
    plot(PCS,cex = 0.75,main = submain,col=rgb(0,0,0,0.5))
  } else {
    if( CAT ){
      plot(PCS,col = factor(DATA[,VAR]),cex = 0.75,
           main = paste0(submain," by ",VAR))
    } else {
      rbPal = colorRampPalette(c("red","blue"))
      plot(PCS,col = rbPal(5)[as.numeric(cut(DATA[,VAR],breaks=5))],
           cex = 0.75,main = paste0(submain," by ",VAR))
    }
  }
}
show_pvalue_hist = function(mat_pvalues,test_type0){
  num_rows = ceiling(sqrt(ncol(mat_pvalues))); #num_rows
  num_cols = ceiling(ncol(mat_pvalues) / num_rows); #num_cols
  par(mfrow=c(num_rows,num_cols),mar=c(5,4,2,1),oma=c(0,0,2,0),bty="n")
  for(cc in seq(ncol(mat_pvalues))){
    hist(mat_pvalues[,cc],xlab="p-value",main=colnames(mat_pvalues)[cc],col="gray")
  }
  mtext(paste0("Linear Model: Type ",test_type0," pvalues"),outer=TRUE,cex=1.4)
  par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,oma=rep(0,4),bty="o")
}
make_dummy = function(x){
  # x = gsub(" ","_",x)
  len_x = length(x)
  all_factors = sort(unique(x))
  num_dummy = length(all_factors) - 1
  fact_matrix = data.frame(matrix(0,nrow=len_x,ncol=max(1,num_dummy)))
  
  if(num_dummy == 0){
    stop("Variable vector is constant")
  } else {
    names(fact_matrix) = paste0(all_factors[-1],"_vs_",all_factors[1])
  }
  
  if(num_dummy == 0){
    for(i in 1:len_x){
      fact_matrix[i,1] = 1
    }
  } else {
    for(i in 1:len_x){
      pos = which(x[i]==all_factors) - 1
      if(is.na(x[i])) fact_matrix[i,] = NA
      else if(pos>0) fact_matrix[i,pos] = 1
    }
  }
  
  fact_matrix
}
smart_ncores = function(){
  ncores = Sys.getenv("SLURM_JOB_CPUS_PER_NODE")
  ncores = ifelse(ncores == "",1,as.numeric(ncores))
  ncores
}
smart_hist = function(x,...){
	hist(x,col="gray",freq=FALSE,...)
	lines(density(x,na.rm=TRUE),lwd=1.5,lty=2,col="blue")
}
bin_cont_var = function(VAR,NUM_GROUPS,ROUND=3,binNUM=FALSE){
  if(FALSE){
    # VAR = sort(runif(50))
    # VAR = aa$AscatM_E
    VAR = aa$LOCI
    VAR = as.numeric(input_mat)
    NUM_GROUPS = 6
    ROUND = 3
    binNUM = TRUE
  }
  
  my_quantiles = as.numeric(quantile(x = VAR,
                                     probs = seq(NUM_GROUPS-1)/NUM_GROUPS,
                                     na.rm = TRUE))
  
  out_VAR = rep(NA,length(VAR))
  for(ii in seq(NUM_GROUPS)){
    if(ii == 1){
      if(binNUM){
        out_VAR[which(VAR <= my_quantiles[ii])] = ii
      } else {
        out_VAR[which(VAR <= my_quantiles[ii])] = paste0(ii,") ",round(min(VAR,na.rm=TRUE),ROUND),
                                                         "-",round(my_quantiles[ii],ROUND))
      }
    } else if(ii == NUM_GROUPS){
      if(binNUM){
        out_VAR[which(VAR > my_quantiles[ii-1])] = ii
      } else {
        out_VAR[which(VAR > my_quantiles[ii-1])] = paste0(ii,") ",round(my_quantiles[ii-1],ROUND),
                                                          "-",round(max(VAR,na.rm=TRUE),ROUND))
      }
    } else {
      if(binNUM){
        out_VAR[which(VAR > my_quantiles[ii-1] & VAR <= my_quantiles[ii])] = ii
      } else {
        out_VAR[which(VAR > my_quantiles[ii-1] & VAR <= my_quantiles[ii])] = paste0(ii,") ",round(my_quantiles[ii-1],ROUND),
                                                                                    "-",round(my_quantiles[ii],ROUND))
      }
    }
  }
  
  # smart_df(VAR,out_VAR)
  if(binNUM) out_VAR = as.character(out_VAR)
  
  out_VAR
}
smart_remove = function(FN,show=FALSE){
	# remove file
	if( file.exists(FN) ){
		blah = file.remove(FN)
	} else if( show ){
		warning(paste0(FN," doesn't exist"))
	}
}
smart_append = function(new_fn,vec_fn){
	# new_fn = name of new file to create
	# vec_fn = vector of files to append together in the order they're stored
	aa = file.create(new_fn)
	bb = sapply(vec_fn,function(x) file.append(new_fn,x))
}
smart_sprintf = function(...){
	orig = sprintf(...)
	orig = gsub("\n","",orig)
	orig = gsub("\t","",orig)
	orig
}



