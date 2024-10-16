# Associations between structural topology of cancer drivers and cancer evolution

Cancer, driven by somatic mutations, disrupts protein function through driver mutations in oncogenes (OCGs) and tumour suppressor genes (TSGs). This study employs structure-based network analysis, a method to understand the structural topology of proteins, to investigate the functional consequences of mutations in cancer driver genes, with a particular focus on lung cancer. By systematically mapping the network properties of each amino acid residue across cancer driver genes using 3D crystallographic proteins structures, and matching these to protein function and cancer mutation frequency in thousands of patients, we found that **structural topology of cancer drivers is strongly linked with cancer evolution. Mutations in highly networked residues were strongly associated with significant functional impairments, highlighting the critical role of central residues in protein stability and function**. In lung cancer, distinct mutation patterns differentiate TSGs, which mutate in highly networked residues in order to lose function, from OCGs, which preferentially mutate in poorly networked residues to retain or gain function, thereby driving cancer. A **novel visualisation method developed using PyMOL functionalities facilitates intuitive mapping of network scores onto protein structures**, revealing regions of high functional importance and providing insights into structural and functional implications of mutations. Future research should aim to expand to all cancer genes. Beyond cancer research, the visualisation tool can adapt to visualise diverse protein metrics and facilitate comparative structural analyses.

Read the [full paper](https://github.com/Yoannli/NetworkAnalysis/blob/main/Final%20Thesis.pdf).

| ![ns](https://github.com/Yoannli/NetworkAnalysis/blob/main/imgs/6GOD_ns.png) |
|:--:| 
| *The wild-type protein of the KRAS gene. On the left, the residue interaction network is shown where nodes represent amino acid residues, edges denote their interactions, and the network score of each residue is colour-coded. On the right, a cartoon depiction of the protein's fold is shown.* |

| ![cbioportal](https://github.com/Yoannli/NetworkAnalysis/blob/main/imgs/cbioportal.png) |
|:--:| 
| *This image represents the culmination of this research. It highlights the association between network scores and mutation frequencies in lung cancer. Tumor-suppressor genes exhibit mutations in highly interconnected residues, whereas oncogenes are found to mutate in less interconnected ones. Additionally, fusion genes and genes with unknown functions are included as control elements.* |


| ![dna](https://github.com/Yoannli/NetworkAnalysis/blob/main/imgs/dna_contact.png) |
|:--:| 
| *Network scores of wild-type p53 protein aligned with DNA contact regions. Residues that are highly networked and functional (173, 176, 274) were located near the DNA contact region. p53 works as a tumour-suppressor gene by binding to DNA and so cancer-causing mutations disrupt its ability to bind DNA.* |
