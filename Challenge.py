import os
from rdflib import Graph, URIRef, Literal, Namespace

working_directory = "/Users/eliascrum/Programs/PythonProjects/"

''' 
Python script for downloading the 'huF7A4DE' VFC file, converting VCF data to RDF format (with splitting for security),
and retrieving 'rs762551 SNP in the CYP1A2 gene' information via a local SPARQL query.

Dependencies:
----
'''


def preprocess_vcf(reference_id):
    """Downloads and Unzips VCF file from online source.

    Args:
        record_id (str): The identity of the file being investigated (sample reference ID)

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
        os.system("gunzip {}".format(reference_id+'.vcf.gz'))
    elif reference_id+'.vcf.gz' in curr_files:
        os.system("gunzip {}".format(reference_id+'.vcf.gz'))

    print("\n\nPreprocessing completed. Continuing to parsing and RDF conversion.")


# basic ontology defined here (just the classes relevant for the challenge - vcf_ontology.ttl)
prefix = Namespace("/Users/eliascrum/Programs/PythonProjects/vcf_ontology.ttl/vcf/Info/")


def parse_to_rdf(rid):
    """Parses VCF file into RDF triples and stores them in different places based on chromosome. Works by
        iterating over each line of an input sample VCF file. Each line contains data for 1 SNP.
        Only nucleotide location (Location), Reference nucleotide (Reference), and Variant nucleotide (Variant) were
        added to RDF data stores (more semantic info could be added but for this challenge it was the only information
        I thought was important.

        Args:
            ?

        Returns:
            VCF sample information in n (where n = number of chromosomes in VCF file) RDF.ttl files
    """
    in_data = open(working_directory + rid + '.vcf', "r").readlines()     # read in vcf data (in chunks maybe?)
    ontology_desc = []
    o = 0
    while in_data[o][0:2] == "##":
        ontology_desc.append(in_data[o])      # for the descriptions of ontology classes ?? Add to ontology doc ??
        o += 1

    prev_chr = 'chr1'       # for dynamic programming and memory efficiency
    curr_chromosome = 1     # to split data into multiple locations
    g = Graph()
    for i in in_data[o+1:]:
        curr_variant = i.split('\t')

        # split data by chromosome for decentralization
        if curr_variant[0] != prev_chr:
            g.serialize(destination=working_directory+"RDF_Data/chromosome_{}.ttl".format(str(curr_chromosome)))
            g = Graph()
            curr_chromosome += 1

        # Important part -- splits each VCF row (SNP) into RDF data and adds to a chromosome-specific graph
        c_id = URIRef("/Users/eliascrum/Programs/PythonProjects/vcf_ontology.ttl/vcf/ID/{}".format(curr_variant[2]))
        g.add((c_id, prefix.Location, Literal(curr_variant[1])))  # ==> position
        g.add((c_id, prefix.Reference, Literal(curr_variant[3])))  # ==> reference base
        g.add((c_id, prefix.Variant, Literal(curr_variant[4])))  # ==> variant base

        # for the last SNP in the VCF file
        if i == in_data[-1]:
            g.serialize(destination=working_directory + "RDF_Data/chromosome_{}.ttl".format(str(curr_chromosome)))

        # tracks current chromosome
        prev_chr = curr_variant[0]

    print("\n\nVCF to RDF Conversion completed. Now beginning SNP information query process.\n")


# global prefix strings to make code look cleaner
info_prefix = "file:///Users/eliascrum/Programs/PythonProjects/vcf_ontology.ttl/vcf/Info/"
id_prefix = "file:///Users/eliascrum/Programs/PythonProjects/vcf_ontology.ttl/vcf/ID/"


def snp_query(snp_id):
    """Parses VCF file into RDF triples and stores them in different places based on chromosome. Works by
        iterating over each line of an input sample VCF file. Each line contains data for 1 SNP.
        Only nucleotide location (Location), Reference nucleotide (Reference), and Variant nucleotide (Variant) were
        added to RDF data stores (more semantic info could be added but for this challenge it was the only information
        I thought was important.

        Args:
            Desired ID of SNP of interest

        Returns:
            ?? Output type ?? about the Reference/Variant information from designated SNP
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
    for kg in os.listdir(working_directory+'RDF_Data'):
        g2 = Graph()
        g2.parse(working_directory+"RDF_Data/" + kg)

        hits = g2.query(q)
        for row in hits:
            for r in row:
                if r is not None:
                    print(r)
                    break
        if done:
            break


# preprocess_vcf("NB72462M")
# parse_to_rdf("NB72462M")
snp_query("rs762551")
