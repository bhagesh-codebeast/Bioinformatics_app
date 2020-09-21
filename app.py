from typing import Counter
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from altair.vegalite.v4.api import Chart
from altair.vegalite.v4.schema.channels import Detail
import streamlit as st
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
from Bio.SeqUtils import MeltingTemp
import matplotlib.pyplot as plt
from stmol import component_3dmol


# APP TITLE
st.title('BIOINFORMATICS TOOLS')

# sidebar contents = https://towardsdatascience.com/creating-streamlit-dashboard-from-scratch-59316a74fa1
st.sidebar.subheader("PICK TOOLS")
select = st.sidebar.selectbox("DROP-DOWN",["ABOUT","DNA FASTA Sequence","Protein FASTA Sequence","Protein PDB File","Protein Data Bank"],key='1')

# OPTION 1 ABOUT
if select == "ABOUT":
    st.header("BIOINFORMATICS SEQUENCE ANALYSIS WEB TOOLS")
    st.subheader("BIOINFORMATICS SEQUENCE ANALYSIS WEB TOOLS")
    st.text("BIOINFORMATICS SEQUENCE ANALYSIS WEB TOOLS")

    #OPTION 2 DNA SEQUENCE
elif select == "DNA FASTA Sequence":

    # main = https://blog.jcharistech.com/2020/05/13/building-a-simple-bioinformatics-app-with-streamlit-and-biopython/
    st.header('DNA SEQUENCE')

    #INPUT USING TEXT BOX
    # seqfileinput = st.text_input('Paste your sequence here',seqfile) = https://blog.jcharistech.com/2019/10/20/streamlit-python-tutorial-crash-course/
    #INPUT USING UPLOAD
    st.set_option('deprecation.showfileUploaderEncoding', False)
    seqfile = st.file_uploader("Upload DNA fasta file",type=["fasta","fa"])

    #BASIC STATS CALCULATION
    if seqfile is not None:
        dnarecord = SeqIO.read(seqfile,"fasta")
        dnaID = dnarecord.id
        dnadescript = dnarecord.description
        dnaseq = dnarecord.seq
        length = len(dnaseq)
        gccont = GC(dnaseq)
        melttemp = MeltingTemp.Tm_GC(dnaseq)
        dnafreq = Counter(dnaseq)
        
        #SETTING RADIO BUTTONS FOR OPTIONS ID|DESCRIPTION|SEQUENCE
        details = st.radio("Sequence Details",("ID","Description","Sequence"))
        if details == "Description":
            st.write(dnadescript)
        elif details == "Sequence":
            st.write(dnaseq)
        elif details == "ID":
            st.write(dnaID)

        #SETTING RADIO BUTTON FOR OPTIONS LENGTH|FREQUENCY TABLE|GC CONTENT|MELTING TEMPERATURE|PLOT NUCLEOTIDE FREQUENCY
        stats = st.radio("Sequence Statistics",("Length","Frequency Table","GC-Content","Melting-Temperature","Plot Nucleotide Frequency"))
        if stats == "Length":
            st.write(length)
        elif stats == "Frequency Table":
            st.write("Frequency Table")
            st.write(pd.DataFrame.from_dict([dnafreq]))
        elif stats == "GC-Content":
            st.write(gccont)
        elif stats == "Melting-Temperature":
            st.write(melttemp)

        #PLOT = https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/examples/sequence_and_translation.html
        elif stats == "Plot Nucleotide Frequency":
            st.write("Pick color")
            adenine = st.beta_color_picker("ADENINE")
            thymine = st.beta_color_picker("THYMINE")
            guanine = st.beta_color_picker("GUANINE")
            cytosine = st.beta_color_picker("CYTOSINE")
            
            #SETTING BUTTON TO PICK COLOUR
            if st.button("Plot"):
                barplot = plt.bar(dnafreq.keys(),dnafreq.values())

                #COLOUR PICKER FOR EACH BASE
                barplot[3].set_color(adenine)
                barplot[2].set_color(thymine)
                barplot[1].set_color(guanine)
                barplot[0].set_color(cytosine)

                #PLOTTING
                st.set_option('deprecation.showPyplotGlobalUse', False)
                st.pyplot()

#OPTION 3 PROTEIN SEQUENCE
elif select == "Protein FASTA Sequence":
    st.header("Protein Sequence")

    #INPUT PROTEIN FASTA FILE
    st.set_option('deprecation.showfileUploaderEncoding', False)
    protfile = st.file_uploader("Upload Protein fasta file",type=["fasta","fa"])

    #CALCULATING BASIC STATS
    if protfile is not None:
        aarecord = SeqIO.read(protfile,"fasta")
        aaID = aarecord.id
        aadescript = aarecord.description
        aaseq = aarecord.seq
        length = len(aaseq)
        aafreq = Counter(aaseq)
        molecw = molecular_weight(aaseq)
        Isoelectric = PA(str(aaseq))
        iso = Isoelectric.isoelectric_point()

        #SETTING BUTTON FOR OPTION PROTEIN ID|PROTEIN DESCRIPTION|PROTEIN SEQUENCE|
        details = st.radio("Sequence Details",("Protein ID","Protein Description","Protein Sequence"))
        if details == "Protein Description":
            st.write(aadescript)
        elif details == "Protein Sequence":
            st.write(aaseq)
        elif details == "Protein ID":
            st.write(aaID)

        #SETTING BUTTON FOR OPTION LENGTH|FREQUENCY TABLE|MOLECULAR WEIGHT|ISO-ELECTRIC POINT|PLOT PROTEIN FREQUENCY
        stats = st.radio("Sequence Statistics",("Length","Frequency Table","molecular weight","Iso-electric Point","Plot Protein Frequency"))
        if stats == "Length":
            st.write(length)
        elif stats == "Frequency Table":
            st.write("Frequency Table")
            st.write(pd.DataFrame.from_dict([aafreq]))
        elif stats == "molecular weight":
            st.write(molecw)
        elif stats == "Iso-electric Point":
            st.write(iso)
        elif stats == "Plot Protein Frequency":
            barplot = plt.bar(aafreq.keys(),aafreq.values())
            st.set_option('deprecation.showPyplotGlobalUse', False)
            st.pyplot()

#OPTION 4 READING PROTEIN PDB FILE
elif select == "Protein PDB File":
    st.write()
    st.write("https://github.com/williamgilpin/pypdb/blob/master/demos/demos.ipynb")
    st.write("https://github.com/williamgilpin/pypdb/blob/master/demos/advanced_demos.ipynb")
    st.write("http://www.wgilpin.com/pypdb_docs/html/")

#OPTION 5 STREAMLIT PDB PLUGIN
elif select == "Protein Data Bank":
    component_3dmol()

#Transcribe or Translate
#1. input 
#2. transcribe
#3. translate
#4. visualise = https://github.com/napoles-uach/streamlit_3dmol

#Alignment
#1. input_number_of_sequence
#2. show aligned sequence
###3. plot alignment = https://github.com/jbkinney/logomaker/tree/master/logomaker/tutorials
