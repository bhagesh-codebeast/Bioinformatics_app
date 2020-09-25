from typing import Counter
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from altair.vegalite.v4.api import Chart
from altair.vegalite.v4.schema.channels import Detail
from google.protobuf import message
import streamlit as st
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
from Bio.SeqUtils import MeltingTemp
import matplotlib.pyplot as plt
from stmol import component_3dmol
from pypdb import *
import requests
import base64


# APP TITLE
st.title('BIOINFORMATICS TOOLS')
st.markdown('******')
# sidebar contents = https://towardsdatascience.com/creating-streamlit-dashboard-from-scratch-59316a74fa1
st.sidebar.subheader("PICK TOOLS")
select = st.sidebar.selectbox("DROP-DOWN",["ABOUT","DNA FASTA Sequence","Protein FASTA Sequence","Protein Data Bank","Protein Visualisation","Table Scraper"],key='1')

# OPTION 1 ABOUT
if select == "ABOUT":
    st.header("LIST OF TOOLS")
    st.subheader("DNA FASTA Sequence")
    st.text("Basic statisttics and visualisation of nucleotide FASTA file")
    st.subheader("Protein FASTA Sequence")
    st.text("Basic statisttics and visualisation of amino acid FASTA file")
    st.subheader("Protein Data Bank")
    st.text("Tool to query and filter Protein Data Bank")
    st.subheader("Protein Visualisation")
    st.text("Tool to visualise proteins using PDB-ID or PDB file")
    st.subheader("Table Scraper")
    st.text("Tool to scrape tables from websites")
    st.markdown('<br><br>', unsafe_allow_html=True)
    st.write(" 👈 Click the arrow on top-left to get started.")

    #OPTION 2 DNA SEQUENCE
elif select == "DNA FASTA Sequence":

    # main = https://blog.jcharistech.com/2020/05/13/building-a-simple-bioinformatics-app-with-streamlit-and-biopython/
    st.header('DNA SEQUENCE')

    #INPUT USING TEXT BOX
    # seqfileinput = st.text_area('Paste your sequence here',"FASTA Sequence") #= https://blog.jcharistech.com/2019/10/20/streamlit-python-tutorial-crash-course/
    # if st.button("Submit"):
    #     seqfile = seqfileinput
    # elif not seqfileinput:
    #     st.warning("Please enter a FASTA Sequence !")
    #     st.stop()

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
elif select == "Protein Data Bank":
    st.write()
    # st.write("https://github.com/williamgilpin/pypdb/blob/master/demos/demos.ipynb")
    # st.write("https://github.com/williamgilpin/pypdb/blob/master/demos/advanced_demos.ipynb")
    # st.write("http://www.wgilpin.com/pypdb_docs/html/")

    #Query PDB Database using pypdb
    st.header("Query Protein Data Bank")
    Options = st.radio("Search By",("pdbid","motif","Author","Litrature","Search By Term","Protein Symmetry","Experimental Method"))
    if Options == "Search By Term":
        searchterm = st.text_input("Enter Search term",'actin network')
        if st.button("Submit"):
            result = Query(searchterm).search()
            st.write(result)
        elif not searchterm:
            st.warning("Please enter a Search term !")
            st.stop()
    elif Options == "pdbid":
        searchterm = st.text_input("Enter PDB ID",'6YYT')
        if st.button("Submit"):
            st.subheader("Protein Description")
            # df = pd.DataFrame.from_dict([describe_pdb(searchterm)])
            # st.write(df.T)
            st.write(describe_pdb(searchterm))
            st.subheader("Protein Information")
            st.write(get_all_info(searchterm))
            # result = get_pdb_file(searchterm, filetype='cif', compression=False)
            # st.write(result)
        elif not searchterm:
            st.warning("Please enter a Search term !")
            st.stop()
    elif Options == "motif":
        searchterm  = st.text_input("Enter Motif",'TAGGY')
        if st.button("Submit"):
            result = Query(searchterm).search()
            st.write(result)
        elif not searchterm:
            st.warning("Please enter a Search term !")
            st.stop()
    elif Options == "Experimental Method":
        searchterm  = st.text_input("Enter experimental method",'SOLID-STATE NMR')
        if st.button("Submit"):
            result = Query(searchterm).search()
            st.write(result)
        elif not searchterm:
            st.warning("Please enter a Search term !")
            st.stop()
    elif Options == "Protein Symmetry":
        searchterm  = st.text_input("Protein Symmetry",'C9')
        if st.button("Submit"):
            result = do_protsym_search(searchterm, min_rmsd=0.0, max_rmsd=1.0)
            st.write(result)
        elif not searchterm:
            st.warning("Please enter a Search term !")
            st.stop()
    elif Options == "Author":
        searchterm  = st.text_input("Author  details",'Perutz, M.F.')
        if st.button("Submit"):
            result = Query(searchterm, query_type='AdvancedAuthorQuery').search()
            st.write(result)
        elif not searchterm:
            st.warning("Please enter a Search term !")
            st.stop()
    elif Options == "Litrature":
        searchterm  = st.text_input("Litrature details",'crispr')
        if st.button("Submit"):
            result = find_papers(searchterm)
            st.write(result)
        elif not searchterm:
            st.warning("Please enter a Search term !")
            st.stop()



#OPTION 5 STREAMLIT PDB PLUGIN
elif select == "Protein Visualisation":
    st.header("Visualise Proteins")
    component_3dmol()

#Table scraper

elif select == "Table Scraper":
    st.header("HTML TABLE SCRAPER")
    try:
        url =  st.text_input("Paste URL Here", value='https://stackexchange.com/leagues/1/alltime/stackoverflow', max_chars=None, key=None, type='default')
        if url:
            arr = ['https://', 'http://']
            if any(crap in url for crap in arr):
                header = {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.169 Safari/537.36",
                "X-Requested-With": "XMLHttpRequest"
                }

                @st.cache(persist=True, show_spinner=False)
                def load_data():
                    html_data = requests.get(url, headers=header)
                    return pd.read_html(html_data.text)

                df = load_data()
                length = len(df)
                    
                if length == 1:
                    st.write("Number of tables: 1 " )
                else: st.write("Number of tables: " + str(length))

                st.subheader("CLICK TO SEE TABLE DETAILS")
                if st.button("Table Details"):
                    st.table(df)
                else: st.empty()

                def createList(r1, r2): 
                    return [item for item in range(r1, r2+1)] 

                st.subheader('SELECT TABLE TO EXPORT')

                r1, r2 = 1, length
                funct = createList(r1, r2)
                ValueSelected = st.selectbox('', funct)
                st.write('You selected table #', ValueSelected)
                df1 = df[ValueSelected -1]


                if df1.empty:
                        st.warning ('This DataFrame is empty!')
                else:
                    df1 = df1.replace(np.nan, 'empty cell', regex=True)
                    st.dataframe(df1)

                    try:
                        csv = df1.to_csv(index=False)
                        b64 = base64.b64encode(csv.encode()).decode()
                        text_display = f'<br><center>Download 👇 Table</center>'
                        st.markdown(text_display, unsafe_allow_html=True)
                        href = f'<center><b><a href="data:file/csv;base64,{b64}" download="filtered_table.csv">Table #{ValueSelected}</a></b></center>'
                        st.markdown(href, unsafe_allow_html=True)
                    except NameError:
                        print('WAIT')
            else:
                st.error ('URL needs to be in a valid format, starting with **https://** or **http://**')
                    
        else:
            st.error("Please Enter URL")
            
    except ValueError:
        st.warning("No Tables Detected")
        
                    

st.markdown("------")
st.markdown("*Made with* 🤘 *by:* [*Bhagesh_Artbeast*](https://github.com/bhagesh-codebeast)")
#Transcribe or Translate
#1. input 
#2. transcribe
#3. translate
#4. visualise = https://github.com/napoles-uach/streamlit_3dmol

#Alignment
#1. input_number_of_sequence
#2. show aligned sequence
###3. plot alignment = https://github.com/jbkinney/logomaker/tree/master/logomaker/tutorials
