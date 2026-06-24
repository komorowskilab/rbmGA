#' @title Permutation Test
#' @author Girish Pulinkala
#' @description Does permutation test to caluclte signifcance of a model
#' @import R.ROSETTA
#' @import data.table
#' @import tibble
#' @import stringr
#' @import grid
#' @import dplyr
#' @import ggplot2
#'
#' @param data decision table used in the training
#' @param iterations A rule set that was generated using the training set
#' @param reducer A character containing name of reducer method: Johnson or Genetic. Default is Johnson.
#' @param discrete Logical. Set TRUE for discrete data. Default is FALSE.
#' @param discreteMethod A character containing discretization method: EqualFrequency, MDL, Naive, SemiNaive or BROrthogonal. Default is EqualFrequency.
#' @param discreteParam A vector containing discretization parameters. May be of different length and values. See examples.
#' @param discreteMask Logical. Set FALSE to disable discretization mask. Default is TRUE.
#' @param cvNum A numeric value of the cross-validation number. Default is 10.
#' @param JohnsonParam A vector containing Johnson reducer parameters.
#' @param GeneticParam A vector containing Genetic reducer parameters. description
#'
#' @examples
#' # example code
#' # df<-fread('../path/to/decision_table)
#' # result<- permutation_test(df, reducer='Johnson',discrete=FALSE)
#' # print(result[[2]])



#Import libraries #FOR IDA
# library(R.ROSETTA)
# library(data.table)
# library(dplyr)
# library(grid)
# library(ggplot2)



permutation_test<- function(data, iterations=100,reducer = 'Johnson',discrete=FALSE,discreteMethod = 'EqualFrequency',
                            cvNum=10,ruleFiltration = TRUE,discreteParam =3,
                            JohnsonParam = list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, Precompute=FALSE, Approximate=TRUE, Fraction=0.95),
                            GeneticParam = list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, Precompute=FALSE, Approximate=TRUE, Fraction=0.95, Algorithm="Simple")){

  shuffle_data<-function(data){

    shuffled_logic <- setNames(list(sample(data[[decision_var]])), decision_var)
    df_shuffled <- do.call(transform, c(list(data), shuffled_logic))

   # df_shuffled=transform( data, data[[decision_var]] = sample(data[[decision_var]]))
    return(df_shuffled)
  }

  col_size<-ncol(data)-1
  decision_var<-names(data)[ncol(data)]
  #sink(file = "lm_output.txt") #should be developed for later
  #message('Orginal labels')
  #message(data[[decision_var]])

  ros<-rosetta(data,reducerDiscernibility = 'Object',discrete = discrete,discreteMethod = discreteMethod ,discreteParam = discreteParam,,reducer = reducer,cvNum = cvNum)
  orginal_accuracy<- ros$quality$accuracyMean
  message('Original classifier built')

  message('Now the shuffling starts')
  quality<-data.frame()

  for(i in 1:iterations){
    set.seed(i) #to generate random initialization and for debugging
    df<-shuffle_data(data)
    ros<-rosetta(df,reducerDiscernibility = 'Object',discrete = discrete,discreteMethod = discreteMethod ,discreteParam = discreteParam,,reducer = reducer,cvNum = cvNum)
    message(paste0('For iteration:',i,' the accuracy is ',ros$quality$accuracyMean))
    quality[i,'Accuracy']<-ros$quality$accuracyMean

  }


  #calculate monte carlo p-values
  r <- sum(ceiling(unlist(quality)*100)/100 >= orginal_accuracy)
  n <- iterations
  # Calculate the adjusted p-value
  p_value <- (r + 1) / (n + 1)

  grob <- grobTree(textGrob(paste0('n=',iterations), x=0.8,  y=0.95, hjust=0,
                            gp=gpar(col="black", fontsize=13, fontface="italic")))

  #confidence intervals
  ci_bounds <- quantile(quality$Accuracy, probs = c(0.025, 0.975), na.rm = TRUE)

  #Max densisty
  dens_x<-max(density(quality$Accuracy)$x)
  dens_y<-max(density(quality$Accuracy)$y)

  permutation_plot<-quality %>%
         ggplot(aes(x = Accuracy)) +
         geom_histogram(aes(y = ..density..), fill = 'royalblue1', color = "white", alpha = 0.8) +
         geom_density(color = "black", linewidth = 1) +
         geom_vline(xintercept = ci_bounds, color = "firebrick", linetype = "dotted", linewidth = 1) +
         geom_vline(aes(xintercept = orginal_accuracy, color = paste0('p<',round(p_value,digits=4))), linewidth = 1.5, linetype = "longdash") +
         scale_x_continuous(breaks = seq(0, 1, 0.1)) +
         theme_classic() +
         annotation_custom(grob) +
         guides(colour = guide_legend(title = 'Original Classifer'))+
         annotate("text",  x = 0.7, y = dens_y * 1.05, label = "Distribution of \nRandom Classifiers", fontface = "italic", color = "darkblue", size = 4) +
         annotate("segment", x = 0.7, y = dens_y*1, xend = dens_x, yend = dens_y * 0.5,  arrow = arrow(length = unit(0.2, "cm")), color = "darkblue")

  print(permutation_plot)

  return(list(quality,permutation_plot))

}




