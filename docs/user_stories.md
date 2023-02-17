Henry is a lab technician with little to no computational experience. He works on a protein design team and needs to develop proteins that are thermostable analogues of existing proteins for industrial processes. He has never used Python or Github, but would benefit from a tool that allows him to generate sequences for wet lab testing. 

Evad is a PI working on computational biology. He has developed highly effective enzymes with little thermostability. He wants a reliable and transparent pipeline for developing high temperature versions of these enzymes. The tool should be well-vetted, open-sourced, have well-understood confidence and limitations, and a clear path for citation.

Rachel has a large library of paired protein sequences but is unsure whether their sequence alignment is indicative of functional similarity. She wants a tool that can quickly parse this dataset to identify which pairs are likely to have similar function.


### Use cases:
1. Sequence functionality matching
Input - SQL database of paired amino acid sequences with species metadata
Output - Score or set of scores mapped to pairwise to input database that is indicative or functional similarity.

2. Thermostable sequence prediction
Input - Mesophilic protein sequence with known (?) functionality.
Output - AA sequence options for thermophilic analogue with similarity scores.

3. Analytic paper trail
Input - Paired AA sequences as in 1
Output - Detailed analysis and data relating to score generation (e.g. 3D folding data, Pfam/GO, etc.)