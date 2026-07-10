set.seed(0)
#Load libraries

if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

if (!require(tibble)) install.packages('tibble')
library(tibble)

if (!require(data.table)) install.packages('data.table')
library(data.table)

if (!require(stringr)) install.packages('stringr')
library(stringr)

if (!require(devtools)) install.packages('devtools')
library(data.table)

if (!require(R.ROSETTA)) install_github('komorowskilab/R.ROSETTA')
library(R.ROSETTA)


combi_gain<-function(data,rules){
  
  col_naming<-function(rules){
    cols<-unique(unlist(strsplit(rules$features, ",")))
    col_names<-vector()
    for(i in cols){
      T<-paste0('P_',i,'_T')
      F<-paste0('P_',i,'_F')
      col_names<-append(col_names,T)
      col_names<-append(col_names,F)
    }
    return(col_names)
  }
  
  names(data)[ncol(data)]<-'decision'
  cuts <- subset(rules,select=c('features','levels','decision',grep('cut', colnames(rules), value=TRUE)))
  colsname<-col_naming(rules)
  df<-data.frame(matrix(nrow=nrow(rules),ncol=length(colsname)))
  names(df)<-colsname
  df <- add_column(df, features<-rules$features, .before = 1)
  df <- add_column(df, levels<-rules$levels, .after = 1)
  df <- add_column(df, decision<-rules$decision, .after=2 )
  names(df)[1:3]<-c('features','levels','decision')
  #no<-melt(table(),value.name = 'no_of_objects')
  no<-melt(data = as.data.table(table(data$decision)),1)
  no$value<-as.numeric(no$value)
  print(no)
  
  ###EXPRS FUnction
  exprs<-function(df,rule_no){
    #print(df)   #checkpoint
    features<-unlist(strsplit(cuts[rule_no,'features'],','))
    #print(features) #checkpoint
    cuts_cond <- unlist(str_split(rules$cuts[rule_no],','))
    #print(cuts_cond) #checkpoint
    cut_counter<-1
    
    for( i in 1:length(features)){
      if(cuts_cond[i] == 'value<cut'){
       # print(features[i]) #checkpoint
        #T CLASS
        n_d_T<-dim(data[data[[features[i]]] < cuts[rule_no,paste0('cut',cut_counter)] & data[['decision']] == cuts[cut_counter,'decision'],])[1] #Number of objects of True decision
        n_u_T<- as.numeric(no[no$V1==cuts[rule_no,'decision'],3])
        #print(cuts[rule_no,'decision'])   #checkpoint
        p_d_T<-n_d_T/n_u_T
        
        print(p_d_T)  
        
        #F CLASS
        n_d_F<-dim(data[data[[features[i]]] < cuts[rule_no,paste0('cut',cut_counter)] & data[['decision']] != cuts[cut_counter,'decision'],])[1] #Number of objects of False decision
        n_u_F<- as.numeric(no[no$V1!=cuts[rule_no,'decision'],3])
        
        p_d_F<-n_d_F/n_u_F
        
        #ASSIGN VALUES
        df[rule_no,paste0('P_',features[i],'_T')]<-p_d_T
        df[rule_no,paste0('P_',features[i],'_F')]<-p_d_F
        
       # print(cuts[rule_no,paste0('cut',cut_counter)])     #checkpoint
        cut_counter<-cut_counter+1
        
      #  print('here')            #checkpoint
      }
      else if(cuts_cond[i] =='cut<value<cut'){
        #T CLASS
        n_d_T<-dim(data[data[[features[i]]] > cuts[rule_no,paste0('cut',cut_counter)] & data[[features[i]]] < cuts[rule_no,paste0('cut',cut_counter+1)] & data[['decision']] == cuts[cut_counter,'decision'],])[1]
        n_u_T<-as.numeric( no[no$V1==cuts[rule_no,'decision'],3])
        p_d_T<-n_d_T/n_u_T
        
        #F CLASS
        n_d_F<-dim(data[data[[features[i]]] > cuts[rule_no,paste0('cut',cut_counter)] & data[[features[i]]] < cuts[rule_no,paste0('cut',cut_counter+1)] & data[['decision']] != cuts[cut_counter,'decision'],])[1]
        print(no[no$V1!=cuts[rule_no,'decision'],3])       
        n_u_F<- as.numeric(no[no$V1!=cuts[rule_no,'decision'],3])
        p_d_F<-n_d_F/n_u_F
        
        #ASSIGN VALUES
        df[rule_no,paste0('P_',features[i],'_T')]<-p_d_T
        df[rule_no,paste0('P_',features[i],'_F')]<-p_d_F
        
       # print(cuts[rule_no,paste0('cut',cut_counter)])     #checkpoint
        cut_counter<-cut_counter+2
        # print('2')        #checkpoint
        
      }
      else if(cuts_cond[i] == 'value>cut'){
        #T CLASS
        n_d_T<-dim(data[data[[features[i]]] > cuts[rule_no,paste0('cut',cut_counter)] & data[['decision']] == cuts[cut_counter,'decision'],])[1]
        n_u_T<- as.numeric(no[no$V1==cuts[rule_no,'decision'],3])
        p_d_T<-n_d_T/n_u_T
        
        #F CLASS
        n_d_F<-dim(data[data[[features[i]]] > cuts[rule_no,paste0('cut',cut_counter)] & data[['decision']] != cuts[cut_counter,'decision'],])[1]
        n_u_F<- as.numeric(no[no$V1!=cuts[rule_no,'decision'],3])
        p_d_F<-n_d_F/n_u_F
        
        #ASSIGN VALUES
        df[rule_no,paste0('P_',features[i],'_T')]<-p_d_T
        df[rule_no,paste0('P_',features[i],'_F')]<-p_d_F
        
       # print(cuts[rule_no,paste0('cut',cut_counter)]). #checkpoint
        cut_counter<-cut_counter+1
        
       # print('3')     #checkpoint
      }
      
      else{
        cat('ERROR')
        break
      }
    }
    return(df)
  }
  
  
  for(k in 1:nrow(rules)){
    df<-exprs(df,k)
    print('1')
  }
  df[is.na(df)] <- 1
  
  th_acc<-function(df,rule_no){
    
    num<-prod(df[rule_no,str_detect(colnames(df),'_T')]) * no[no$V1==df[1,'decision'],3]
    den<-( prod(df[rule_no,str_detect(colnames(df),'_T')]) *  no[no$V1==df[1,'decision'],3] ) + (prod(df[rule_no,str_detect(colnames(df),'_F')]) * no[no$V1!=df[1,'decision'],3])
    th_acc <- num/den
    df[rule_no,'Therotical_Accuracy']<-th_acc
    
    return(df)
  }
  
  for(k in 1:nrow(rules)){
    df<-th_acc(df,k)
  }
  
  #rules_gain<-rules %>% left_join(dplyr::select(df,features,levels,Therotical_Accuracy,grep('cut',colnames(rules))),by=c('features'='features','levels'='levels')) %>% dplyr::select(c(features,levels,decision,supportRHS,accuracyRHS,Therotical_Accuracy,grep('cut',colnames(rules))))
  rules_gain<-rules %>% left_join(dplyr::select(df,features,levels,Therotical_Accuracy),by=c('features'='features','levels'='levels')) %>% dplyr::select(c(features,levels,decision,supportRHS,accuracyRHS,Therotical_Accuracy,grep('cut',colnames(rules))))
  rules_gain$Combinatorial_gain <- (rules_gain$Therotical_Accuracy-rules_gain$accuracyRHS)
  return(rules_gain)
  
}


#EXAMPLE 
#PLEASE USE BINARY DECISION, YET to be developed for multiple
#ros<-rosetta(autcon)
#rules<-ros$main
#rules_gain<-combi_gain(data= autcon,rules)
#vis<-visunet(rules_gain[rules_gain$Therotical_Accuracy > 0.8,])


# 
# data<-fread("../data/Rosetta_Decision_Tables/D7.csv")
# data<-data[data$protectionStatus=='Prot' | data$protectionStatus=='NonProt',]
# data<-data[data$Group =='S' | data$Group =='O' | data$Group =='ABL1_a',]
# ros<-readRDS('../data/Rosetta_results/gen.rds')
# rules<-ros$main



#rules<-readRDS('../data/Rosetta_results/MPC_classifier.rds')
# rules_gain<-combi_gain(data,classifier)
# rules_gain<-distinct(rules_gain)
#fwrite(rules_gain,'../data/Rosetta_results/rules_w_combigain.csv')

#rules_df<-rules_gain[rules_gain$Therotical_Accuracy > 0.8 & rules_gain$Combinatorial_gain > -0.2,]

#vis<-visunet(rules_gain[rules_gain$Therotical_Accuracy > 0.8,])
# vis<-visunet(rules_gain[rules_gain$Combinatorial_gain>-0.2,])



