import os
import argparse
from pathlib import Path
from rdflib import Graph, URIRef, Literal, Namespace


working_directory = Path('/Users/eliascrum/Programs/PythonProjects')

run = argparse.ArgumentParser(description='Download, parse, and/or query VCF files.')
run.add_argument('-s', '--sample', nargs='?', default='NB72462M',
                 help='-d flag specifies a Personal Genome Project sample ID you would like to download ' +
                      'and/or parse to RDF')
run.add_argument('-q', '--query', nargs='?', default='rs762551',
                 help='-q flag specifies a SNP ID and will return information about it')
run.add_argument('-t', '--test', action='store_true',
                 help='runs a quick test to ensure proper instillation')
args = run.parse_args()

# sample = vars(args)['sample']
# query = vars(args)['query']
test = True #vars(args)['test']
sample = 'NB72462M'
query = 'rs762551'
''' 
Python script for downloading the 'huF7A4DE' VFC file, converting VCF data to RDF format (with splitting for security),
and retrieving 'rs762551 SNP' information via a local SPARQL query.

Dependencies:
    Linux OS terminal (this script throws many errors in windows command prompt)
    Argparse - Python Package
    RDFLib - Python Package
    
----
'''


def preprocess_vcf(reference_id):
    """Downloads and Unzips VCF file from online source.

    Args:
        reference_id (str): The identity of the file being investigated (sample reference ID)

    Returns:
        A downloaded, unzipped VCF file: a .vcf file that can be parsed
    """

    curr_files = os.listdir(working_directory)
    if (reference_id+'.vcf' or reference_id+'.vcf.gz') not in curr_files:
        os.system(
            "wget --mirror --no-parent --no-host --cut-dirs=1 " +
            "https://531155966bc06bca5de62439c00ce64b-282.collections.ac2it.arvadosapi.com/_/" +
            "{}".format(reference_id+'.vcf.gz')
        )
        # Linux systems
        os.system("gunzip {}".format(reference_id+'.vcf.gz'))

    elif reference_id+'.vcf.gz' in curr_files:
        os.system("gunzip {}".format(reference_id+'.vcf.gz'))

    print("\n\n\tPreprocessing completed. Continuing to parsing and RDF conversion.")


# basic ontology defined here (just the classes relevant for the challenge - vcf_ontology.ttl)
prefix = Namespace(working_directory / "vcf_ontology.ttl" / "vcf" / "Info")


def parse_to_rdf(rid):
    """Parses VCF file into RDF triples and stores them in different places based on chromosome. Works by
        iterating over each line of an input sample VCF file. Each line contains data for 1 SNP.
        Only nucleotide location (Location), Reference nucleotide (Reference), and Variant nucleotide (Variant) were
        added to RDF data stores (more semantic info could be added but for this challenge it was the only information
        I thought was important.

        Args:
            Sample VCF file ID

        Returns:
            VCF sample information in n# of RDF.ttl files (where n = number of chromosomes in VCF file)
    """
    in_data = open(working_directory / (rid + '.vcf'), "r").readlines()     # read in vcf data (in chunks maybe?)
    ontology_desc = []
    o = 0
    while in_data[o][0:2] == "##":
        ontology_desc.append(in_data[o])      # for the descriptions of ontology classes ?? Add to ontology doc ??
        o += 1

    prev_chr = 'chr1'       # for dynamic programming and memory efficiency
    curr_chromosome = 'chr1'     # to split data into multiple locations
    g = Graph()

    # iterates over each line of the VCF file and converts data to RDF
    for i in in_data[o+1:]:
        curr_variant = i.split('\t')

        # because we are only working with known SNPs for this challenge (so ignoring those without SNP IDs)
        if curr_variant[2] != '.':

            # split data by chromosome for decentralization
            if curr_variant[0] != prev_chr:
                g.serialize(destination=working_directory / "RDF_Data" / "{}.ttl".format(str(curr_chromosome)))
                g = Graph()
                curr_chromosome = curr_variant[0]

            # Important part -- splits each VCF row (SNP) into RDF data and adds to a chromosome-specific graph
            c_id = URIRef(str(working_directory / "vcf_ontology.ttl" / "vcf" / "ID" / "{}".format(curr_variant[2])))
            g.add((c_id, prefix.Location, Literal(curr_variant[1])))  # ==> position
            g.add((c_id, prefix.Reference, Literal(curr_variant[3])))  # ==> reference base
            g.add((c_id, prefix.Variant, Literal(curr_variant[4])))  # ==> variant base

            # for the last SNP in the VCF file
            if i == in_data[-1]:
                g.serialize(destination=working_directory / "RDF_Data" / "{}.ttl".format(str(curr_chromosome)))

            # tracks current chromosome
            prev_chr = curr_variant[0]

    print("\n\tVCF to RDF Conversion completed. Now beginning SNP information query process.\n")


# global prefix strings to make code look cleaner -- LINUX (need Windows version here) *****
info_prefix = "file://{}".format(working_directory / "vcf_ontology.ttl" / "vcf" / "Info")
id_prefix = "file://{}".format(working_directory / "vcf_ontology.ttl" / "vcf" / "ID") + '/'



def snp_query(snp_id):
    """Parses VCF file into RDF triples and stores them in different places based on chromosome. Works by
        iterating over each line of an input sample VCF file. Each line contains data for 1 SNP.
        Only nucleotide location (Location), Reference nucleotide (Reference), and Variant nucleotide (Variant) were
        added to RDF data stores (more semantic info could be added but for this challenge it was the only information
        I thought was important.

        Args:
            Desired ID of SNP of interest

        Returns:
            SNP_Information.txt file that contains the Variant allele information from designated SNP
    """

    # SPARQL Query to retrieve designated SNP variant allele nucleotide
    q = """
    PREFIX inf: <%s>
    PREFIX id: <%s>
    SELECT ?p ?o
    WHERE {
        id:%s inf:Variant ?o .
    }""" % (info_prefix, id_prefix, snp_id)

    # perform SPARQL query on each chromosome's RDF graph until we get the result
    done = False
    for kg in os.listdir(working_directory / 'RDF_Data'):
        g2 = Graph()
        g2.parse(working_directory / "RDF_Data" / kg)

        hits = g2.query(q)
        output = open(working_directory / 'SNP_Information.txt', 'w')
        for row in hits:
            for r in row:
                if r is not None:
                    output.write('The SNP %s is found on %s and contains the variant "%s" allele'
                                 % (snp_id, kg[:-4], str(r)))
                    done = True
                    break
        if done:
            output.close()
            print('\nSuccess! The SNP information you are looking for can be found in the "SNP_Information.txt" file.\n')
            break


# test case
if test:
    print("You selected to run the test version. If there are any errors they will be displayed below.")
    preprocess_vcf('test')
    parse_to_rdf('test')
    snp_query('rs62637819')

else:
    preprocess_vcf(sample)
    parse_to_rdf(sample)
    snp_query(query)
