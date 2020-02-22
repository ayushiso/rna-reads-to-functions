import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax

up_file = sys.argv[1]
down_file = sys.argv[2]
onto_type = sys.argv[3]

ontology = pd.read_csv("go_slim_mapping.tab", sep='\t', header=None)
ontology = ontology[[0, 3, 4]]
ontology.columns = ['gene', 'o_type', 'onto']

# combine ontology and ge set and get mean log2fold change values
def get_onto_groups(gene_file):
    genes = pd.read_csv(gene_file, dtype=str, sep='\t')
    genes = genes.set_index('gene')
    gene_ont = pd.merge(genes, ontology, how='inner', on='gene')
    gene_ont = gene_ont.loc[gene_ont['o_type'] == onto_type]
    gene_ont["log2FoldChange"] = pd.to_numeric(gene_ont["log2FoldChange"])
    l2fc = gene_ont.groupby('onto')['log2FoldChange'].mean()

    # get gene ontology categories with highest representation
    gene_counts = gene_ont["onto"].value_counts()
    category_counts = ontology["onto"].value_counts()
    proportions = pd.merge(gene_counts, category_counts, how='inner', left_index=True, right_index=True)
    proportions['prop'] = proportions['onto_x']/proportions['onto_y']
    proportions = pd.merge(proportions, l2fc,how='inner', left_index=True, right_index=True)
    proportions = proportions.sort_values('prop', ascending=False)[:10]

    return proportions

# plotting code
up_proportions = get_onto_groups(up_file)
down_proportions = get_onto_groups(down_file)
proportions = pd.concat([up_proportions, down_proportions])

plt.figure(0, figsize=(10,11))
plt.grid(True)
plt.scatter(proportions.index.values, proportions['log2FoldChange'], s = proportions['prop']*5000)
plt.xticks(rotation=90)
vars = up_file.split(".")[0].split("/")[1].split("_")
if len(vars) == 3:
    plt.title("Differentially expressed genes with " + vars[1]  + ": " + vars[2], fontsize=20)
else:
    plt.title("Differentially expressed genes with " + vars[1], fontsize=20)
plt.xlabel("Gene ontology", fontsize=20)
plt.ylabel("Log2fc value", fontsize=20)
plt.tight_layout()
if len(vars) == 3:
    plt.savefig(vars[1] + "_" + vars[2] +  "_" + onto_type + ".png")
else:
    plt.savefig(vars[-1] + "_" + onto_type + ".png")
