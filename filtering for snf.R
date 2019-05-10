
load("/Users/tesssenftle/Downloads/methyldata.Rda")
#write.csv(transposedmeth, "/Users/tesssenftle/Downloads/dnameth.csv")
install.packages("SNFtool")
library(SNFtool)
#transposed = t(mRNA)
#write.csv(transposed, "/Users/tesssenftle/Downloads/mrna.csv")
#write.csv(dnameth_filtered_na_filtered, "/Users/tesssenftle/Downloads/dnamethfilteredfinal.csv")

transposedcn <- read.csv(file = "/Users/tesssenftle/Downloads/copynumber.csv")
transposedcn2 <- transposedcn[,-1]
rownames(transposedcn2) <- transposedcn[,1]

transposed = read.csv(file="/Users/tesssenftle/Downloads/mrna.csv")
transposedmeth = read.csv(file="/Users/tesssenftle/Downloads/dnameth.csv")
transposedmeth2 <- transposedmeth[,-1]
rownames(transposedmeth2) <- transposedmeth[,1]

transposed2 <- transposed[,-1]
rownames(transposed2) <- transposed[,1]



cropped1 = dmeth2[0:200]
cropped2 = fcommongenes[0:200]
cropped3 = copynumber_filtered[0:200,0:200]

cropped4 = dnameth_filtered[0:200,0:200]
cropped5 = finalcopynumber_filtered[0:200,0:200]
cropped6 = finalde[0:200]


cropped8 = finalmrna_filtered

cropped9 = finalcopynumber_filtered[0:414,0:500]

tempmrna <- transposed2
tempdnameth <- transposedmeth2
tempcopynumber <- transposedcn2

tempmrna2 <- tempmrna
tempdnameth2 <-tempdnameth 

copynumberp = rownames(transposedcn2)
dnamethp = rownames(transposed2)
mrnap = rownames(transposedmeth2)

dnap = rownames(tempdnameth2)
rnap = rownames(tempmrna2)

 
copynumberp2 = substr(mrnap, 1, 12)
mrnap2= substr(mrnap, 1, 12)
dnamethp2 = substr(dnamethp, 1, 12)

mrnap3= substr(dnap, 1, 12)
dnamethp3 = substr(rnap, 1, 12)

head(copynumberp2)
head(mrnap2)
head(dnamethp2)
common = intersect(dnamethp2,mrnap2)
#length(common)

#common = intersect(dnamethp3,mrnap3)


#temprna <-unique(tempmrna)
#tempdna <-unique(tempdnameth)





 uniquecopynumber <- unique(copynumberp2)
 uniquemrna <- unique(mrnap2)
 uniquednameth <- unique(dnamethp2)
# 
# 
# 
# tempmrna <- transposed2
# tempdnameth <- transposedmeth2
# 
# tempmrnanames <- rownames(tempmrna)
# tempmrnanames2 <- substr(tempmrnanames, 1, 12)
# 
# uniquetempmrna <-unique(tempmrnanames2)
# 
# rownames(tempmrna) <- uniquetempmrna
# 
# tempdnanames <- rownames(tempdnameth)
# tempdnanames2 = substr(tempdnanames, 1, 12)
# 
# uniquetempdna <-unique(tempdnanames2)

.rowNamesDF(tempmrna, make.names=TRUE) <- uniquemrna
.rowNamesDF(tempdnameth, make.names=TRUE) <- uniquednameth
.rowNamesDF(tempcopynumber, make.names=TRUE) <- uniquecopynumber
write.csv(tempcopynumber, "/Users/tesssenftle/Downloads/tempcopynumber.csv") 


tempcopynumber <- read.csv(file = "/Users/tesssenftle/Downloads/tempcopynumber.csv")
#next step
#finalcopynumber <- tempcopynumber[,-1]
#rownames(finalcopynumber) <- tempcopynumber[,1]

mrna_filtered = subset(as.data.frame(finalmrna), rownames(finalmrna) %in% common)
dnameth_filtered = subset(as.data.frame(tempdnameth), rownames(tempdnameth) %in% common)
copynumber_filtered = subset(as.data.frame(finalcopynumber), rownames(finalcopynumber) %in% common)

#write.csv(copynumber_filtered,"/Users/tesssenftle/Downloads/copynumber_filtered.csv")
#write.csv(mrna_filtered, "/Users/tesssenftle/Downloads/mrna_filtered.csv") 
#write.csv(dnameth_filtered, "/Users/tesssenftle/Downloads/dnameth_filtered.csv") 

dnameth_filtered <- read.csv(file="/Users/tesssenftle/Downloads/dnamethfilteredfinal.csv")

copynumber_filtered <- read.csv(file="/Users/tesssenftle/Downloads/copynumber_filtered.csv")
mrna_filtered <- read.csv(file="/Users/tesssenftle/Downloads/mrna_filtered.csv")

finaldnameth_filtered  <- dnameth_filtered [,-1]
rownames(finaldnameth_filtered) <- dnameth_filtered[,1]
write.csv(finaldnameth_filtered, "/Users/tesssenftle/Downloads/finaldnameth_filtered.csv")


finalcopynumber_filtered2  <- finalcopynumber_filtered [,-1]
rownames(finalcopynumber_filtered2) <- finalcopynumber_filtered[,1]
write.csv(finalcopynumber_filtered, "/Users/tesssenftle/Downloads/finalcopynumber_filtered.csv")

finalmrna_filtered  <- mrna_filtered [,-1]
rownames(finalmrna_filtered) <- mrna_filtered[,1]
write.csv(finalmrna_filtered, "/Users/tesssenftle/Downloads/finalmrna_filtered.csv")
testmrna_filtered <- read.csv(file="/Users/tesssenftle/Downloads/finalmrna_filtered.csv")

#ncropped1 <- cropped1[order(row.names(cropped1)),] 
#ncropped2 <- cropped2[order(row.names(cropped2)),] 

## First, set all the parameters:
K = 10;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 10; 	# Number of Iterations, usually (10~20)

## Data1 is of size n x d_1, where n is the number of patients, d_1 is the number of genes, e.g.
## Data2 is of size n x d_2, where n is the number of patients, d_2 is the number of methylation, e.g.
#data(Data1)
#data(Data2)

## Here, the simulation data (Data1, Data2) has two data types. They are complementary to each other. And two data types have the same number of points. The first half data belongs to the first cluster; the rest belongs to the second cluster.

#truelabel = c(matrix(1,100,1),matrix(2,100,1)); ##the ground truth of the simulated data;


## Calculate distance matrices(here we calculate Euclidean Distance, you can use other distance, e.g,correlation)

## If the data are all continuous values, we recommend the users to perform standard normalization before using SNF, though it is optional depending on the data the users want to use.  

Data1 = standardNormalization(cropped7);
Data2 = standardNormalization(cropped8);
Data3 = standardNormalization(cropped9);


## Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows; if the data is discrete, we recommend the users to use ""
Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
Dist3 = dist2(as.matrix(Data3),as.matrix(Data3));

## next, construct similarity graphs
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)
W3 = affinityMatrix(Dist3, K, alpha)

## These similarity graphs have complementary information about clusters.
#displayClusters(W1,truelabel);
#displayClusters(W2,truelabel);

## next, we fuse all the graphs
## then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF(list(W1,W2, W3), K, T)

## With this unified graph W of size n x n, you can do either spectral clustering or Kernel NMF. If you need help with further clustering, please let us know. 
## for example, spectral clustering

C = 5			# number of clusters
group = spectralClustering(W, C); 	# the final subtypes information

## you can evaluate the goodness of the obtained clustering results by calculate Normalized mutual information (NMI): if NMI is close to 1, it indicates that the obtained clustering is very close to the "true" cluster information; if NMI is close to 0, it indicates the obtained clustering is not similar to the "true" cluster information.

displayClusters(W, group);
#SNFNMI = calNMI(group, truelabel)

## you can also find the concordance between each individual network and the fused network

ConcordanceMatrix = concordanceNetworkNMI(list(W, W1,W2,W3),C);


#dnameth_filtered_na_filtered <- dnameth_filtered[,which(!is.na(apply(dnameth_filtered, 2, mean)))]



cluster_outcome <- as.data.frame(group, row.names(finalcopynumber_filtered))











################################################################################
# We also provide an example using label propagation to predict the labels of new data points below.
# How to use SNF with multiple views

# Load views into list "dataL"
# load("Digits.RData")
data(Digits)
# Set the other parameters
K = 20 # number of neighbours
alpha = 0.5 # hyperparameter in affinityMatrix
T = 20 # number of iterations of SNF
# Normalize the features in each of the views (optional)
# dataL = lapply(dataL, standardNormalization)

# Calculate the distances for each view
distL = lapply(dataL, function(x) dist2(x, x))

# Construct the similarity graphs
affinityL = lapply(distL, function(x) affinityMatrix(x, K, alpha))
################################################################################
# An example of how to use concordanceNetworkNMI

Concordance_matrix = concordanceNetworkNMI(affinityL, 3);

## The output, Concordance_matrix,  shows the concordance between the fused network and each individual network. 

################################################################################
# Example of how to use SNF to perform subtyping
# Construct the fused network
W = SNF(affinityL, K, T)
# perform clustering on the fused network.
clustering = spectralClustering(W,3);
# use NMI to measure the goodness of the obtained labels.
NMI = calNMI(clustering, label);

################################################################################
# Provide an example of predicting the new labels with label propagation

# Load views into list "dataL" and the cluster assignment into vector "label"
data(Digits)

# Create the training and test data
n = floor(0.8*length(label)) # number of training cases
trainSample = sample.int(length(label), n)
train = lapply(dataL, function(x) x[trainSample, ]) # Use the first 150 samples for training
test = lapply(dataL, function(x) x[-trainSample, ]) # Test the rest of the data set
groups = label[trainSample]

# Set the other
K = 20
alpha = 0.5
t = 20
method = TRUE

# Apply the prediction function to the data
newLabel = groupPredict(train,test,groups,K,alpha,t,method)

# Compare the prediction accuracy
accuracy = sum(label[-trainSample] == newLabel[-c(1:n)])/(length(label) - n)

