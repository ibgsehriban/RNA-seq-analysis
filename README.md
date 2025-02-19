#Step1: Create Environment for analysis

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA
import os

#Step2: Load RNA-seq count data
os.chdir("/Users/sehriban/Desktop/d63")  # Change  working directory
df = pd.read_csv("counts.tsv", sep="\t") ## you have to change the file name to your file name
print(df.shape)
print(df.head())

#Step3: Data Preprocessing
#Filter lowly expressed genes
min_count = 10  # Minimum count threshold
df = df.loc[df.sum(axis=1) > min_count]

#Select samples and genes of interest
new_df= df.iloc[:, [6, 7, 3, 4]]  # Python uses 0-based indexing
df1 = new_df.dropna()  # Removes any rows with NA values

# Log-transform data for visualization
df2 = np.log1p(df1)

#Step5: Create metadata
metadata = pd.DataFrame({
    "Sample": ["SRR6506142.bam", "SRR6506143.bam", "SRR6506139.bam", "SRR6506140.bam"],
    "Condition": ["Control", "Control", "Treatment", "Treatment"]
})

# Save metadata as a CSV file
metadata.to_csv("metadata.csv", index=False)

# Display metadata
print(metadata)

#Step6: Create Plots
# Visualize expression distribution
# Boxplot of log-transformed RNA-seq counts
plt.figure(figsize=(10, 6)) # Set figure size
sns.boxplot(data=df2) # Create boxplot
plt.xticks(rotation=90) # Rotate x-axis labels
plt.title("Log-transformed RNA-seq Count Distribution") # Add title
plt.show() # Show plot

# Normalize data using TPM
def calculate_tpm(counts): # First create function to calculate TPM
    gene_lengths = counts.sum(axis=1)  # Assume sum of counts represents length
    rpk = counts.div(gene_lengths, axis=0)
    tpm = rpk.div(rpk.sum()) * 1e6  # Normalize to 1 million reads per sample        
    return tpm # TPM values for each gene

data_tpm = calculate_tpm(df2)  # Calculate TPM
data_tpm = data_tpm.dropna()  # Remove rows with NaN values
print(data_tpm.head()) # Display first 5 rows


# Step7: Perform PCA on the data
pca = PCA(n_components=2) # Initialize PCA
pca_result = pca.fit_transform(data_tpm.T)  
pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"], index=new_df.columns)        
pca_df = pca_df.merge(metadata, left_index=True, right_on="Sample", how="left")


print(pca_df.head())  # Ensure "PC1" and "PC2" exist

plt.figure(figsize=(8, 6))

# Ensure Condition is correctly assigned
print(pca_df.head())  

sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Condition")

plt.title("PCA of RNA-seq Data")
plt.show()


#Step8: Differential Expression Analysis
condition1 = metadata[metadata["Condition"] == "Control"]["Sample"] # Control samples
condition2 = metadata[metadata["Condition"] == "Treatment"]["Sample"]   # Treatment samples

t_stat, p_values = ttest_ind(df1[condition1], df1[condition2], axis=1)
diff_expr_results = pd.DataFrame({"Gene": df1.index, "p_value": p_values})
diff_expr_results["-log10(p_value)"] = -np.log10(diff_expr_results["p_value"])

mean_control = df1[condition1].mean(axis=1)
mean_treatment = df1[condition2].mean(axis=1)
log_fc = np.log2(mean_treatment + 1) - np.log2(mean_control + 1)
diff_expr_results["log2_FC"] = log_fc

#Volcano Plot of Differential Expression
plt.figure(figsize=(8, 6))
sns.scatterplot(data=diff_expr_results, x=diff_expr_results.index, y="-log10(p_value)")
plt.axhline(y=-np.log10(0.05), color='red', linestyle='--')
plt.title("Volcano Plot of Differential Expression")        # Add title
plt.xlabel("Genes")       # Add x-axis label
plt.ylabel("-log10(p-value)")           # Add y-axis label
plt.show()


