#clean the environment
rm(list=ls(all=TRUE))
nrep = 500
#read in the datasets
setkeep = data.frame(read_csv('/Users/macszerez.com/Desktop/portfolio/conor/data/setkeep.csv'))
reskeep = data.frame(read_csv('/Users/macszerez.com/Desktop/portfolio/conor/data/reskeep.csv'))

# first we clean reskeep
reskeep=na.omit(reskeep)
reskeep

#set up counters, p_values, estimates and standard errors
counter=reskeep[,36:38]
#p_value sequence
ip=c(seq(3,27,3))
#estimate sequence
ie=c(seq(1,25,3))
#standard error sequence
is=c(seq(2,26,3))
# adjusted r2 sequence
ir = 28:35
ir
reskeep_p_values=cbind(reskeep[,ip],counter)
reskeep_p_values
reskeep_est=cbind(reskeep[,ie],counter)
reskeep_est
reskeep_se=cbind(reskeep[,is],counter)
reskeep_se
reskeep_R2=cbind(reskeep[,ir],counter)
reskeep_R2

#----------------------------------------
#-----------setkeep cleaning-------------
#----------------------------------------

setkeep=na.omit(setkeep[,1:10])
setkeep=as.data.frame(setkeep)
ii = length(setkeep$par_a)
setkeep=cbind(setkeep,NA)
colnames(setkeep)=c("nmz","ndz","par_a",'par_c', 'par_e', 'par_g', 'par_b', 'par_x', 'p_pgs', 'p_A','set')
#it gives a set id to each set so when we join it to the result, it will be clear which set beyonds to which results
for(i in 1:ii){
  setkeep[i,11]=i
}
# sequence setkeep to be able to join it to the results
setkeep_full <- setkeep[rep(seq_len(nrow(setkeep)), each = nrep), ]
#check if all good
head(reskeep)
check=as.data.frame(setkeep_full[,11] == reskeep$set)
colnames(check) = 'BOOL'
check_BOOL=check[which(check$BOOL == FALSE),]
#should return an empty list
if (length(check_BOOL) == 0) {
  print('The set IDs match')
} else {
  print(paste(length(check_BOOL), 'problems in the dataset'))
}
#merging setkeep and reskeep for final dataframes

final_df_p = cbind(setkeep_full,reskeep_p_values)
final_df_est = cbind(setkeep_full, reskeep_est)
final_df_se = cbind(setkeep_full,reskeep_se)
#the set is included twice due to joining, we delete that
final_df=cbind(setkeep_full,reskeep)[,-11]

#-------------------------------------------------------
#--------------------POWER ANALYSIS---------------------
#-------------------------------------------------------


#number of rows
nrowdf=nrow(final_df_p)
nrowdf/rep
# Split the original dataset into smaller datasets
row_indices <- seq_len(nrow(final_df_p))
nrep <- ceiling(nrow(final_df_p) / rep)
df_list <- split(final_df_p, rep(1:nrep, each = rep, length.out = nrow(final_df_p)))
dim(df_list$`1`)
# Calculate the power for each column in each dataset
alpha <- 0.05
power_list <- lapply(df_list, function(x) {
  apply(x, 2, function(y) {
    sum(y < alpha) / length(y)
  })
})
# merge the data together
power_table <- do.call(rbind, power_list)
power_table = power_table[,12:(ncol(power_table)-3)]
power_table
#merge the data with setkeep
final_power_table=cbind.data.frame(setkeep,power_table)
final_power_table=round(final_power_table,8)
#save the power table
final_power_table
#write.csv(final_power_table,"/Users/macszerez.com/Desktop/VU GBH/Genes in Health and Behaviour 1st year/Internship I/results/set1/final_power_table_set1.csv")

#calculate the mean power for each test among all the values
means_set=apply(final_power_table,2,mean)
means_set
means_set=as.data.frame(t(means_set))
means_set=means_set[,12:20]
transformed_set=t(means_set)
rownames(transformed_set)=colnames(means_set)
colnames(transformed_set)='power'
transformed_set
#save the mean power for each test
#write.csv(transformed_set,"/Users/macszerez.com/Desktop/VU GBH/Genes in Health and Behaviour 1st year/Internship I/results/set2/mean_of_power_set2.csv")

# get the average adjusted R-square values for the tests using the same techniques as with the power analysis
df_list_R <- split(reskeep_R2, rep(1:nrep, each = rep, length.out = nrow(reskeep_R2)))
R2_list=lapply(df_list_R,colMeans)
R2_list
R2_table=do.call(rbind, R2_list)[,1:8]
R2_table
final_R2_table=cbind.data.frame(setkeep,R2_table)
final_R2_table
#write.csv(final_R2_table,"/Users/macszerez.com/Desktop/VU GBH/Genes in Health and Behaviour 1st year/Internship I/results/set1/R2_table_set1.csv")
mean_overall=as.data.frame(apply(final_R2_table,2,mean))
mean_overall=t(t(mean_overall[12:20,]))
mean_overall
rownames(mean_overall)=colnames(final_R2_table)
mean_overall

final_power_table

# subset the power result to a strictly dizigotyc dataset
power_values_dz=final_power_table[13:16]
# and a monozygotic-dizygotic dataset
power_values_dzmz=final_power_table[17:20]

#df_se=data.frame(read.csv('/Users/macszerez.com/Desktop/VU GBH/Genes in Health and Behaviour 1st year/set_1_full.csv'))
#df_se

#Visualization of the results

library(ggplot2)
library(ggpubr)

# Define the parameter values (x-axis)
parameter_values <- 1:15

# Define colors for the power values
power_colors <- c("red", "blue", "green", "orange")

# Create a data frame with parameter, power values, and colors
data <- data.frame(Parameter = rep(parameter_values, each = 4),
                   Power = as.vector(t(power_values_dz)),
                   Models = rep(1:4, times = 15))

# Create a scatter plot with colored points
dz=ggplot(data, aes(x = Parameter, y = Power, color = factor(Models))) +
  geom_line(size = 1.5) +
  scale_color_manual(values = power_colors,
                     labels = c('M0dz','M1dz','M2dz','M3dz')) +
  scale_x_continuous(breaks = parameter_values, labels = parameter_values) +
  xlab("Parameter") +
  ylab("Power") +
  ggtitle("Change in Power with Parameter Values (set 3 DZ)")


# Create a data frame with parameter, power values, and colors
data <- data.frame(Parameter = rep(parameter_values, each = 4),
                   Power = as.vector(t(power_values_dzmz)),
                   Models = rep(1:4, times = 15))

# Create a scatter plot with colored points
mz=ggplot(data, aes(x = Parameter, y = Power, color = factor(Models))) +
  geom_line(size = 1.5) +
  scale_color_manual(values = power_colors,
                     labels = c('M0dzmz','M1dzmz','M2dzmdz','M3dzmz')) +
  scale_x_continuous(breaks = parameter_values, labels = parameter_values) +
  xlab("Parameter") +
  ylab("Power") +
  ggtitle("Change in Power with Parameter Values (set 3 DZ-MZ)")
mz
dz
figure <- ggarrange(dz,mz,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
figure

# Check the average power per model in DZ and MZ-DZ dataset

barchart <- data.frame(Category = rownames(transformed_set),
                       Value = as.vector(t(transformed_set[, 1])))

ggplot(barchart, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("red", "blue", "green", "orange", "purple",
                               "cyan", "magenta", "yellow", "brown"))
geom_hline(yintercept = 0.75, color = "blue", linetype = "dashed", size = 1)
                    
  
  

