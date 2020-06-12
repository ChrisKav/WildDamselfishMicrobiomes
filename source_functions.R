plot_ordered_bar <- function (physeq, x = "Sample",
                              y = "Abundance",
                              fill = NULL,
                              leg_size = 0.5,
                              title = NULL) {
  require(ggplot2)
  require(phyloseq)
  require(plyr)
  require(grid)
  bb <- psmelt(physeq)
  samp_names <- aggregate(bb$Abundance, by=list(bb$Sample), FUN=sum)[,1]
  .e <- environment()
  bb[,fill]<- factor(bb[,fill], rev(sort(unique(bb[,fill])))) #fill to genus
  bb<- bb[order(bb[,fill]),] # genus to fill
  p = ggplot(bb, aes_string(x = x, y = y,
                            fill = fill),
             environment = .e, ordered = FALSE)
  p = p +geom_bar(stat = "identity",
                  position = "stack",
                  color = "black")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  p = p + guides(fill = guide_legend(override.aes = list(colour = NULL), reverse=TRUE)) +
    theme(legend.key = element_rect(colour = "black"))
  p = p + theme(legend.key.size = unit(leg_size, "cm"))
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt,
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}
se <- function(x) sqrt(var(x)/length(x))
summarize_taxa = function(physeq, Rank, GroupBy = NULL){
  Rank <- Rank[1]
  if(!Rank %in% rank_names(physeq)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){
    GroupBy <- GroupBy[1]
    if(!GroupBy %in% sample_variables(physeq)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt
  mdt = fast_melt(physeq)
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    sdt = data.table(SampleID = sample_names(physeq),
                     var1 = get_variable(physeq, GroupBy))
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  Nsamples = nsamples(physeq)
  summarydt = mdt[, list(meanRA = sum(RelativeAbundance)/Nsamples,
                         sdRA = sd(RelativeAbundance),
                         seRA = sd(RelativeAbundance)/sqrt(nsamples(physeq)),
                         seRA = se(RelativeAbundance),
                         minRA = min(RelativeAbundance),
                         maxRA = max(RelativeAbundance)),
                  by = c(Rank, GroupBy)]
  return(summarydt)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' @title outputs a FASTA file from a phyloseq object
#'
#' @description This function outputs a FASTA-formatted text file from a \code{phyloseq} object
#'
#' @param ps A \code{phyloseq} object that contains \code{\link[phyloseq]{refseq}}.
#' If there the \code{refseq} slot is not filled, this function will try pull the
#' sequences from \code{\link[phyloseq]{get_taxa}}
#'
#' @param file (optional) A file name that ends in ".fasta" or ".fa".
#' If a file name is not supplied, the file will be named after the phyloseq object.
#'
#' @param rank (optional) A taxonomic rank from the \code{\link[phyloseq]{tax_table}} which will be used to name the sequences.
#' If no rank is supplied, samples will be named \code{ASV_#}
#'
#' @return This function saves a FASTA-formatted text file from the input \code{phyloseq} object.
#'
#' @examples
#' save_fasta(ps)
#' save_fasta(ps = ps, file = "sequences.fasta", rank = "Genus")

save_fasta <- function(ps = ps, file = NULL, rank = NULL){
  
  if(is.null(ps)){
    message("Phyloseq object not found.")
  }
  
  if(is.null(file)){
    file <- paste0(deparse(substitute(ps)), ".fasta")
  }
  
  if(is.null(rank) | !rank %in% rank_names(ps)){
    message("Rank not found. Naming sequences sequentially (i.e. ASV_#).")
    seq_names <- paste0("ASV_", 1:ntaxa(ps))
  } else {
    seq_names <- make.unique(unname(tax_table(ps)[,rank]), sep = "_")
  }
  
  if(!is.null(refseq(ps))){
    seqs <-as.vector(refseq(ps))
  } else{
    message("refseq() not found. Using taxa names for sequences.")
    if(sum(grepl("[^ACTG]", rownames(tax_table(ps)))) > 0){
      message("Taxa names do not appear to be DNA sequences. Proceed with caution.")
    }
    seqs <-get_taxa(ps)
  }
  
  for (i in 1:ntaxa(ps)){
    cat(paste(">", seq_names[i], sep=""), file=file, sep="\n", append=TRUE)
    cat(seqs[i], file=file, sep="\n", append=TRUE)
  }
  message(paste0(ntaxa(ps), " sequences written to <", file, ">."))
}

#' Calculate Good's Coverage
#'
#' Calculates Good's coverage from a community data matrix with samples as rows and OTUs as columns.
#'
#' @param com a vegan compatible community data matrix.
#'
#' @return A table with the headings number of singletons, number of sequences, and Good's coverage for each sample in rows.
#' @export
#'
#' @references Good, I. J. 1953. The Population Frequencies of Species and the Estimation of Population Parameters. Biometrika 40:237-264.
#'
#' @examples
#'
goods <-
function(com){
  no.seqs <- rowSums(com)
  sing <- com==1
  no.sing <- apply(sing, 1, sum)
  goods <- 100*(1-no.sing/no.seqs)
  goods.sum <- cbind(no.sing, no.seqs, goods)
  goods.sum <- as.data.frame(goods.sum)
  return(goods.sum)
}

save.excel <-function(.list, default = 'var1', path = ''){
    require("XLConnect")
    .name <- as.list(match.call())[2]
    if(is.language(.name[[1]])) wb_name <- paste0(paste0(path, default, collapse = '/'), '.xlsx')
    if(is.symbol(.name[[1]])) wb_name <- paste0(paste0(path, as.character(.name), collapse = '/'), '.xlsx')
    wb <- loadWorkbook(wb_name, create = TRUE)
    createSheet(wb, names(.list))
    writeWorksheet(wb,.list, names(.list),header=FALSE)
    saveWorkbook(wb)
    }

core.fit=function(otu_table) {
  otu_table_no_tax=otu_table[-length(otu_table)]
  otu_table_no_tax=t(otu_table_no_tax)
  ##Calculate the average of total individuals per community source
  
  N <- mean(apply(otu_table_no_tax, 1, sum))
  
  ##Calculate the average relative abundance of each taxa across communities
  p.m <- apply(otu_table_no_tax, 2, mean)
  p.m <- p.m[p.m != 0]
  p <- p.m/N
  
  ##Calculate the occurrence frequency of each taxa across communities
  freq <- apply(otu_table_no_tax>0, 2, mean)
  freq <- freq[freq != 0]
  
  
  ##Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  ##Calculate the limit of detection  
  d = 1/N
  
  ##Goodness of fit for Poisson model
  require(Hmisc)
  pois.pred <- ppois(d, N*p, lower.tail=F)
  pois.pred.ci<- binconf(pois.pred*nrow(otu_table), nrow(otu_table), alpha=0.95, method="exact", return.df=TRUE)
  
  otu_matrix=cbind(freq,p,pois.pred,pois.pred.ci)
  
  SS_res=sum((freq*log(freq/pois.pred))-(freq-pois.pred))
  SS_tot=sum((freq*log(freq/mean(freq))-(freq-mean(freq))))
  
  Rsqr.pois<- 1 - (SS_res/SS_tot)
  
  # R square  adj
  SS_res2=sum((freq*log(freq/pois.pred))-(freq-pois.pred))+(length(subset(otu_matrix,freq>=pois.pred.ci[,1]))/2)
  
  SS_tot2=sum((freq*log(freq/mean(freq))-(freq-mean(freq))))
  
  
  Rsqr.pois.adj <- 1 - (SS_res2/SS_tot2)
  
  # Root Mean Squared Error
  RMSE.pois =sum((pois.pred - freq)^2)/(length(freq) - 1)
  
  #p-value
  t=(mean(freq)-pois.pred)/(sd(freq)/sqrt(length(freq)))
  p_value=2*pt(-abs(t),df=length(freq)-1)
  
  
  ## Split the the microbial core table and the variable community table
  core.table_Poisson=subset(otu_table,row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=pois.pred.ci[,2])))
  core.table_30=subset(otu_table,row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.3)))
  core.table_40=subset(otu_table,row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.4)))
  core.table_50=subset(otu_table,row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.5)))
  core.table_60=subset(otu_table,row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.6)))
  core.table_70=subset(otu_table,row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.7)))
  core.table_80=subset(otu_table,row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.8)))
  core.table_90=subset(otu_table,row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.9)))
  core.table_100=subset(otu_table,row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=1.0)))
  
  
  nocore.table_Poisson=subset(otu_table,!row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=pois.pred.ci[,2])))
  nocore.table_30=subset(otu_table,!row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.3)))
  nocore.table_40=subset(otu_table,!row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.4)))
  nocore.table_50=subset(otu_table,!row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.5)))
  nocore.table_60=subset(otu_table,!row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.6)))
  nocore.table_70=subset(otu_table,!row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.7)))
  nocore.table_80=subset(otu_table,!row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.8)))
  nocore.table_90=subset(otu_table,!row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=0.9)))
  nocore.table_100=subset(otu_table,!row.names(otu_table) %in% row.names(subset(otu_matrix,freq>=1.0)))
  
  return(list(core_Poisson=core.table_Poisson,core_30=core.table_30,core_40=core.table_40,core_50=core.table_50,core_60=core.table_60,core_70=core.table_70,core_80=core.table_80,core_90=core.table_90,core_100=core.table_100,nocore_Poisson=nocore.table_Poisson,nocore_30=nocore.table_30,nocore_40=nocore.table_40,nocore_50=nocore.table_50,nocore_60=nocore.table_60,nocore_70=nocore.table_70,nocore_80=nocore.table_80,nocore_90=nocore.table_90,nocore_100=nocore.table_100,otu_matrix=otu_matrix,Rsqr.pois.adj=Rsqr.pois.adj,p_value=p_value))
}
