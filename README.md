
Data_visualiszation_Omic shinyapp


📘 How to Use the Shiny RNA-seq / Omics Explorer
This app lets you explore RNA-seq or other omics data interactively, creating publication-ready heatmaps, PCA plots, and gene expression boxplots.
access the app:
https://janangawra.shinyapps.io/Data_visualiszation_Omic/



1️⃣ Uploading Files
Counts / normalized expression file: Upload a CSV, TSV, TXT, or Excel file (genes × samples).
Gene ID column: Select the column containing gene identifiers.
Sample metadata (colData): Upload metadata with one row per sample (conditions, replicates, etc.).
Sample ID column: Select the column that uniquely identifies each sample (must match column names in the count matrix).
2️⃣ Coloring Options
Select a metadata factor to color plots.
Manually pick colors for each group.
Optionally highlight a single sample.
3️⃣ Heatmap Tab
Choose all genes, top variable genes, or a custom list.
Select scaling: none, Z-score rows, or columns.
Enable/disable clustering for rows/columns.
Click 'Render heatmap' to visualize.
4️⃣ PCA Tab
Performs PCA on log2-transformed and optionally scaled data.
Colors samples by the selected metadata factor.
Edit axis labels and titles.
5️⃣ Gene Expression Tab
Type or select a gene ID to view expression.
Choose plot type: Boxplot, Violin, or Jitter.
Log-transform and show points optionally.
📂 Example Input Format
Counts file:
Gene_ID, CD_CTRL_1, CD_CTRL_2, CD_C1_1, CD_C1_2
ENSDARG0000000018, 123, 245, 312, 221
ENSDARG0000000321, 432, 534, 312, 298
...

colData file:
SampleName, Condition, Replicate
CD_CTRL_1, Control, R1
CD_C1_1, Low_dose, R1
CD_C3_1, Medium_dose, R1
CD_C5_1, High_dose, R1
