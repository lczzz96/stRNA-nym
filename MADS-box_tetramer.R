library(dplyr)
library(ggsci)
library(ggplot2)
library(cowplot)
library(patchwork)
library(plotrix)
library(ggforce)
library(ggpubr)
library(RColorBrewer)
library(plotrix)
library(pheatmap)
library(grid)
library(reshape2)

is_valid = function(w, x, y, z) {
  w >= 0 && x >= 0 && y >= 0 && z >= 0 &&           
    2 * w + x <= ap1 &&                             
    x + y <= ap3 &&                                
    x + y <= pi &&                                
    y + 2 * z <= ag &&                             
    2 * w + x + y + 2 * z <= sep                    
}

optimal_flower_parts2 = function(ap1, ap3, pi, ag, sep) {
  max_total = 0
  best_combinations = list()
  
  for (w in 0:floor(ap1/2)) {  
    for (x in 0:(ap1 - 2*w)) { 
      for (y in 0:min(ap3 - x, pi - x)) {  
        max_z_ag = floor((ag - y) / 2)       
        max_z_sep = floor((sep - 2*w - x - y) / 2) 
        max_z = min(max_z_ag, max_z_sep)
        for (z in 0:max_z) {
          if (is_valid(w, x, y, z)) {
            total = w + x + y + z
            if (total > max_total) {
              max_total = total
              best_combinations = list(list(w = w, x = x, y = y, z = z, total = total))
            } else if (total == max_total) {
              best_combinations = c(best_combinations, list(list(w = w, x = x, y = y, z = z, total = total)))
            }
          }
        }
      }
    }
  }
  
  if (length(best_combinations) > 0) {
    cat("Optimal combination(s) found:\n")
    for (sol in best_combinations) {
      cat(sprintf("Sepal = %d, Petal = %d, Stamen = %d, Carpel = %d (Total = %d)\n", 
                  sol$w, sol$x, sol$y, sol$z, sol$total))
    }
    
    if (length(best_combinations) > 1) {
      avg_w = mean(sapply(best_combinations, function(sol) sol$w))
      avg_x = mean(sapply(best_combinations, function(sol) sol$x))
      avg_y = mean(sapply(best_combinations, function(sol) sol$y))
      avg_z = mean(sapply(best_combinations, function(sol) sol$z))
      avg_total = mean(sapply(best_combinations, function(sol) sol$total))
      
      cat("\nAverage of all optimal solutions:\n")
      cat(sprintf("Sepal = %.2f, Petal = %.2f, Stamen = %.2f, Carpel = %.2f (Total = %.2f)\n", 
                  avg_w, avg_x, avg_y, avg_z, avg_total))
      
      cat("\nDetailed statistics:\n")
      cat(sprintf("Number of optimal solutions: %d\n", length(best_combinations)))
      cat(sprintf("Range of Sepal: %d - %d\n", 
                  min(sapply(best_combinations, function(sol) sol$w)),
                  max(sapply(best_combinations, function(sol) sol$w))))
      cat(sprintf("Range of Petal: %d - %d\n", 
                  min(sapply(best_combinations, function(sol) sol$x)),
                  max(sapply(best_combinations, function(sol) sol$x))))
      cat(sprintf("Range of Stamen: %d - %d\n", 
                  min(sapply(best_combinations, function(sol) sol$y)),
                  max(sapply(best_combinations, function(sol) sol$y))))
      cat(sprintf("Range of Carpel: %d - %d\n", 
                  min(sapply(best_combinations, function(sol) sol$z)),
                  max(sapply(best_combinations, function(sol) sol$z))))
    }
  } else {
    cat("No valid combinations found.\n")
  }
}

### each stamen vs IT vs OT ####
mads_table=read.table("../data/mads_box.table",row.names = 1, header=T)

for(i in 1:nrow(mads_table)) {
  mads_table[i,]=as.integer(mads_table[i,]*100)
}

for(i in 1:nrow(mads_table)) {
  print(i)
  ap1=mads_table[i,1]
  ap3=mads_table[i,2]
  pi=mads_table[i,3]
  ag=mads_table[i,4]
  sep=mads_table[i,5]
  
  print(c(ap1, ap3, pi, ag, sep))
  best_combinations=optimal_flower_parts2(ap1, ap3, pi, ag, sep)
}

prop_matrix = prop_matrix %>%
  mutate_all(as.numeric)
rownames(prop_matrix) = paste0("Sample", 1:6)

prop_matrix_long = prop_matrix %>%
  mutate(Sample = rownames(.)) %>%
  melt(id.vars = "Sample", variable.name = "Variable", value.name = "Value")

prop_matrix_long$Sample = factor(prop_matrix_long$Sample, levels = rev(c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")))

ggplot(prop_matrix_long, aes(x = Variable, y = Sample, fill = Value)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = round(Value, 3)), color = "black", size = 3) +
  scale_fill_gradient2(low = "white", high = "red", 
                       midpoint = max(prop_matrix_long$Value)/2,
                       name = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))


ggplot(prop_matrix_long, aes(x = Variable, y = Sample, fill = Value)) +
  geom_tile(color = "darkgray", linewidth = 0.2) +  # 改为黑色边线
  scale_fill_gradientn(colors = c("white", "orange2", "tomato"),
                       values = scales::rescale(c(0, 0.2, 0.5, 1)),
                       limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme_void() +
  theme(legend.position = "none")

### all stamens vs IT vs OT ####
stamen_counts = colSums(mads_table[1:4,])
mads_table2 = rbind(mads_table, stamen = stamen_counts)
mads_table2 = mads_table2[c(7,5,6),]

for(i in 1:nrow(mads_table2)) {
  mads_table2[i,]=as.integer(mads_table2[i,]*100)
}

for(i in 1:nrow(mads_table2)) {
  print(i)
  ap1=mads_table2[i,1]
  ap3=mads_table2[i,2]
  pi=mads_table2[i,3]
  ag=mads_table2[i,4]
  sep=mads_table2[i,5]
  
  print(c(ap1, ap3, pi, ag, sep))
  best_combinations=optimal_flower_parts2(ap1, ap3, pi, ag, sep)
}


### all stamens vs IT vs OT (bulk data)####
mads_table3 = read.table("../data/mads_box_bulk.table", sep = "\t", header = T, row.names = 1)

for(i in 1:nrow(mads_table3)) {
  mads_table3[i,]=as.integer(mads_table3[i,]*10)
}

for(i in 1:nrow(mads_table3)) {
  print(i)
  ap1=mads_table3[i,1]
  ap3=mads_table3[i,2]
  pi=mads_table3[i,3]
  ag=mads_table3[i,4]
  sep=mads_table3[i,5]
  
  print(c(ap1, ap3, pi, ag, sep))
  best_combinations=optimal_flower_parts2(ap1, ap3, pi, ag, sep)
}
