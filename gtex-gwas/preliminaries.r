
##library(tidyverse) 


## Connect to DB
##install.packages('RPostgreSQL')
library('RPostgreSQL')
drv <- dbDriver("PostgreSQL")
db.gwas.v6p  <- dbConnect(drv, host="gene2pheno.ccsriudvs1y0.us-east-1.rds.amazonaws.com", port="5432", dbname="public_gwas_v6p_hapmap_2017_12_04", user="metaxcan_ro", password="M3t4Xc4nR0")
db.uk.v6p <- dbConnect(drv, host="gene2pheno.ccsriudvs1y0.us-east-1.rds.amazonaws.com", port="5432", dbname="ukb_neale_2017_10_14", user="metaxcan_ro", password="M3t4Xc4nR0")


## Query functions
query0 = paste0(  
  "SELECT "
  ," g.gene_name"
  ,", m.zscore"
  ,", m.effect_size"
  ,", m.pval"
  ,", p.tag as phenotype"
  ,", t.tissue as tissue"
  ,", m.pred_perf_R2"
  ,", m.pred_perf_pval"
  ,", m.pred_perf_qval"
  ,", m.n_snps_used"
  ,", m.n_snps_model"
  ,", mi.p_smr"
  ,", mi.p_heidi"
  ,", mi.coloc_p0"
  ,", mi.coloc_p1"
  ,", mi.coloc_p2"
  ,", mi.coloc_p3"
  ,", mi.coloc_p4"	
  ,", tr.p as twas_pvalue "
  ,", g.gene"
  ," FROM gene AS g "
  ," INNER JOIN metaxcan_result AS m ON g.id = m.gene_id "
  ," INNER JOIN tissue AS t ON t.id = m.tissue_id "
  ," INNER JOIN pheno AS p ON p.id = m.pheno_id "
  ," LEFT JOIN metaxcan_result_info as mi on m.id = mi.metaxcan_result_id "
  ," LEFT JOIN twas_relationship as trr ON m.id = trr.metaxcan_result_id "
  ," LEFT JOIN twas_result as tr ON trr.twas_result_id = tr.id "    
)

query.pheno = function(pheno,Rthres=0.,fdb=db)
{
  query =  paste0(query0, " WHERE p.tag like '" , pheno  ,"%' " )
  if(Rthres>0) query = paste0(query, " and m.pred_perf_R2 >= " , Rthres)
  dbGetQuery(fdb,query) %>% arrange(pval)
}

query.tissue = function(tissue,Rthres=0.,fdb=db)
{
  query =  paste0(query0," WHERE t.tissue like '" , tissue,  "%' " )
  if(Rthres>0) query = paste0(query, " and m.pred_perf_R2 >= ", Rthres)
  dbGetQuery(fdb,query)
}

query.gene = function(genename,Rthres=0.,fdb=db)
{
  query =  paste0(query0," WHERE g.gene_name like '", genename,  "%' " )
  if(Rthres>0) query = paste0(query, " and m.pred_perf_R2 >= ", Rthres)
  dbGetQuery(fdb,query) %>% arrange(pval)
}

#
query_multi_tissue_base <- paste0( "SELECT ",
                                   " p.tag as phenotype,",
                                   " g.gene,",
                                   " g.gene_name,",
                                   " mt.pvalue,",
                                   " mt.n,",
                                   " mt.n_indep,",
                                   " mt.p_i_best,",
                                   " mt.p_i_worst,",
                                   " mt.eigen_max,",
                                   " mt.eigen_min",
                                   " FROM gene AS g ",
                                   " INNER JOIN multi_tissue_result AS mt ON g.id = mt.gene_id ",
                                   " INNER JOIN pheno AS p ON p.id = mt.pheno_id " )

query.pheno.mult = function(pheno,fdb=db)
{
  query =  paste0(query_multi_tissue_base, " WHERE p.tag like '" , pheno  ,"%' " )
  dbGetQuery(fdb,query) %>% arrange(pvalue)
}

query.gene.mult = function(genename,Rthres=0.,fdb=db)
{
  query =   paste0(query_multi_tissue_base," WHERE g.gene_name like '", genename,  "%' " )
  dbGetQuery(fdb,query) %>% arrange(pvalue)
}

### --------

## drop repeated vars but keeping first row for given var
droprep = function(df,repvar,sortvar,decreasing=FALSE)
{
  if(repvar %in% names(df))
  {  
    if(sortvar %in% names(df))
    {
      tab = table(df[[repvar]])
      nametab = names(tab)
      uniquelist = nametab[tab==1]
      replist =  nametab[tab>1]
      ind = df[[repvar]] %in% uniquelist
      res = df[ind,]
      tempo = df[!ind,]
      tempo = tempo[order(tempo[[sortvar]],decreasing=decreasing),]
      for(vv in replist)
      {
        rowa = tempo[tempo[[repvar]]==vv,][1,]
        res = rbind(res,rowa)
      }
      return(res)
    } else {warning(sortvar %&% ' (sortvar) not in data frame')}
  } else {warning(repvar %&% ' (repvar) not in data frame')}
}




