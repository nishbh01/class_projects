import streamlit as st
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter
import io
import pandas as pd
import seaborn as sns

# Data viz
import matplotlib
matplotlib.use("Agg")

def codon_frequency(sequence, codons):
    frequency = {}
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        
        if codon in codons:
            if codon in frequency:
                frequency[codon] += 1
            else:
                frequency[codon] = 1
    return frequency


isoleucine_codons = ['ATA', 'ATC', 'ATT']
arginine_codons = ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]
alanine_codons = ["GCA", "GCC", "GCG", "GCT"]
asparagine_codons = ["AAC", "AAT"]
aspartic_acid_codons = ["GAC", "GAT"]
cysteine_codons = ["TGC", "TGT"]
glutamic_acid_codons = ["GAA", "GAG"]
glutamine_codons = ["CAA", "CAG"]
glycine_codons = ["GGC", "GGT", "GGA", "GGG"]
histidine_codons = ["CAC", "CAT"]
leucine_codons = ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"]
lysine_codons = ["AAA", "AAG"]
methionine_codons = ["ATG"]
phenylalanine_codons = ["TTC", "TTT"]
proline_codons = ["CCC", "CCT", "CCA", "CCG"]
serine_codons = ["TCC", "TCT", "TCA", "TCG", "AGC", "AGT"]
threonine_codons = ["ACC", "ACT", "ACA", "ACG"]
tryptophan_codons = ["TGG"]
tyrosine_codons = ["TAC", "TAT"]
valine_codons = ["GTA", "GTC", "GTG", "GTT"]

codon_lists = {
    'Isoleucine': isoleucine_codons,
    'Asparagine': asparagine_codons,
    'Glutamic acid': glutamic_acid_codons,
    'Alanine': alanine_codons,
    'Arginine': arginine_codons,
    'Aspartic acid': aspartic_acid_codons,
    'Cysteine': cysteine_codons,
    'Glutamine': glutamine_codons,
    'Glycine': glycine_codons,
    'Histidine': histidine_codons,
    'Leucine': leucine_codons,
    'Lysine': lysine_codons,
    'Methionine': methionine_codons,
    'Phenylalanine': phenylalanine_codons,
    'Proline': proline_codons,
    'Serine': serine_codons,
    'Threonine': threonine_codons,
    'Tryptophan': tryptophan_codons,
    'Tyrosine': tyrosine_codons,
    'Valine': valine_codons
}

# Main App Coding in Streamlit (st)

def main():
	""" A simple bioinformatics app"""
	
	st.title("Relative Codon Usage")
	menu = ["Intro", "Codon Analysis", "Gene Clustering", "About"]
	choice = st.sidebar.selectbox("Select Activity", menu)
	if choice == "Intro":
		st.write(
                         """
        Welcome to my WEB APP 'Codoverse' 

        While translation of transcribed bases into the amino acids is a straightforward phenomenon in multicellular organisms, 
        it is found in recent studies that translation might be affected by the type of codon used for a specific amino acid. 
        There are 64 different codons, and 61 of them code for functional codons – other than stop codons. 
        Only 20 amino acids are coded by these 64 different amino acids, leading to an environment of higher fidelity and variability. 
        It was long thought that synonymous codons for an amino acid are used indiscriminately, but there have been frequent studies suggesting otherwise.

        **This study focused on:**
        - Analyzing the usage pattern of synonymous codons, especially isoleucine codons, in breast cancer transcripts compared to stable control genes (APOE4).
        - (But you will be able to upload fasta files for any genes)
        - Studying mutations in cancer-related genes, including promoters and suppressors, to identify statistically significant transitions in nucleotide preferences for isoleucine codons.
        - Analyzing tRNAs from cancer cells to determine if there is a preference for a specific codon among the isoleucine codons.

        A web-based app will be developed to accompany this research, allowing researchers to upload their transcript files and examine codon preferences.
        The uploaded FASTA file can contain one or multiple sequences.

        ### Get Started
        To get started, select an activity from the sidebar menu.
        """
    )
	elif choice == "Gene Clustering":
         st.subheader("Generate Clustermaps Based On Codon Usage")
         uploaded_files = st.file_uploader("Upload FASTA File", type=["fasta", "fa"], accept_multiple_files=True)
         if uploaded_files:
            codon_menu = ['Isoleucine', 'Asparagine', 'Glutamic acid', 'Alanine', 'Arginine', 'Aspartic acid', 'Cysteine', 'Glutamine', 'Glycine', 'Histidine', 'Leucine', 'Lysine', 'Methionine', 'Phenylalanine', 'Proline', 'Serine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine']
            sel_codon = st.selectbox("Which codon's usage are you interested in?", codon_menu, key='codon_selection')
            freq_dfs = {}
            for uploaded_file in uploaded_files:
                    gene_name = uploaded_file.name
                    gene_name = gene_name.replace('.fasta', '')
                    try:
                        fasta_sequences = [str(record.seq) for record in SeqIO.parse(io.StringIO(uploaded_file.getvalue().decode('utf-8')), "fasta")]
                    except:
                        fasta_sequences = SeqIO.read(uploaded_file, "fasta")
                        fasta_sequences = fasta_sequences.seq
                    frequency_cal = codon_frequency(str(fasta_sequences), codon_lists[sel_codon])
                    freq_dfs[gene_name] = frequency_cal

            def cluster_maps():
                        df = pd.DataFrame.from_dict(freq_dfs, orient="index")
                        combined_df = pd.DataFrame(freq_dfs).T.reset_index().rename(columns={'index':'Gene'})
                        combined_df = pd.melt(combined_df, id_vars=['Gene'], var_name='synonymous codon', value_name='count')
                        combined_df = combined_df.fillna(0)

                        # Create row and column label DataFrames
                        pivot_df = combined_df.pivot(index='Gene', columns='synonymous codon', values='count')
                        st.write(df)
                        # Create the cluster map
                        sns.clustermap(pivot_df, row_cluster=True, col_cluster=True, cmap='coolwarm', linewidths=.5, figsize=(10,10))
                        plt.title('Relative Evolutionary Distances Between Genes Based on codon usage')
                        st.pyplot()

            if st.button("Cluster Genes"):
                cluster_maps()  # calls function when user clicks.
	elif choice == "Codon Analysis":
		st.subheader("Bias in the Use of Synonymous codons")
		input_file = st.file_uploader("Upload FASTA File", type=["fasta", "fa"])
		if input_file is not None:
			try:
				sequences1 = [str(record.seq) for record in SeqIO.parse(io.StringIO(input_file.read().decode('utf-8')), "fasta")]
			except:
				sequences1 = SeqIO.read(input_file, "fasta")
				sequences1 = sequences1.seq
                
			codon_menu = ['Isoleucine', 'Asparagine', 'Glutamic acid', 'Alanine', 'Arginine', 'Aspartic acid', 'Cysteine', 'Glutamine', 'Glycine', 'Histidine', 'Leucine', 'Lysine', 'Methionine', 'Phenylalanine', 'Proline', 'Serine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine']

			def calculate_frequency():
				codon_list = codon_lists[sel_codon]
				frequency = codon_frequency(str(sequences1), codon_list)
				st.write(sel_codon)
				st.write((f"{sel_codon}: {frequency}"))
				codons = frequency.keys()
				count = frequency.values()
				plt.bar(codons, count, color = ['r', 'b', 'g', 'c'])
				plt.show()
				st.set_option('deprecation.showPyplotGlobalUse', False)
				st.pyplot() # shows the visualization in web.
            
			sel_codon = st.selectbox("Which codon's usage are you interested in?", codon_menu, key='codon_selection')
			st.button("Calculate Codon Frequency", on_click= calculate_frequency) # calls function when user clicks.
			
			
	elif choice == "About":
              st.image("/Users/nischal/Library/CloudStorage/GoogleDrive-nbhanda1@ramapo.edu/My Drive/undergrad/spring 2023/AdvBio/project/bioinfoapp/IMG_1053.jpg", caption="Nischal Bhandari", use_column_width=True)
    
              st.write(
        """
        ## About Me

       *Hello!*
        My name is Nischal Bhandari, and I am a bioinformatics major at Ramapo College of New Jersey. Currently, I am working on my senior project, which focuses on studying Codon Usage in Breast Cancers. I am passionate about understanding the molecular mechanisms underlying cancer development and progression. 

        ## Curiosity and Future Plans

        My curiosity lies in exploring the complex interplay between genetics, epigenetics, and environmental factors in cancer. I am particularly interested in deciphering the role of specific genetic mutations and their impact on cancerous processes. By integrating bioinformatics tools and computational analysis, I aim to uncover novel insights that could contribute to the development of more effective diagnostic and therapeutic strategies.

        """
    )
              st.markdown("[Download Resume](https://drive.google.com/file/d/1xZCXVDJcvxMeF-79zCPJWowoLaH9hGD6/view?usp=sharing)")
              # Embed the PDF file using an iframe
              import requests
              st.markdown(f'<iframe src="{"https://drive.google.com/file/d/1xZCXVDJcvxMeF-79zCPJWowoLaH9hGD6"}" width="100%" height="600"></iframe>', unsafe_allow_html=True)


if __name__ == "__main__":
	main()

