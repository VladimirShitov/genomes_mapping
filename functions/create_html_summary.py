import dominate
from dominate import tags
from dominate.util import raw  # For escaping html in strings

from constants import ORGANISM


def create_html_summary(n_genomes, n_clusters, n_big_clusters):
    doc = dominate.document(title='Genomes mapping summary')

    with doc.head:
        tags.style("""
        body {
            background-color: #e1f5fe;
            font-family: sans-serif;
            font-size: 18px;
        }
        .container {
            width: 80%;
            margin: 0 auto;
            background-color: white;
            padding: 15px 30px;
        }
        span {
            font-weight: 600;
        }
        cite {
            border-left: solid black 3px;
            font-size: 12px;
        }
        .big_img {
            width: 100%;
        }
        """)

    with doc.body:
        with tags.div(cls='container'):
            tags.h1('Creating the total list of genes for {}'.format(ORGANISM))

            tags.p(raw('Total list has <span>{} genomes</span>'.format(n_genomes)
                       + ' and <span>{} clusters</span>.'.format(n_clusters)))
            tags.img(src='../plots/clusters_sizes.png', alt='clusters sizes')

            tags.p(raw('There are <span>{} clusters</span>'.format(n_big_clusters) + ' of <span>size > 5</span>.'))
            tags.img(src='../plots/n_clusters_by_threshold.png', alt='N clusters by threshold')

            tags.p(raw('Check the <i>results/</i> folder for files, which you can use for the following analysis:'))

            with tags.ul():
                tags.li(raw('<i><b>proteome.faa</b></i> — fasta file with all the protein coding genes from'
                            + 'all the genomes. Gene names are concatenated with their genome ID. There is also'
                            + 'some extra information in brackets. Example of the header:<br> '
                            + '<cite>>WP_000002298.1[genome:GCF_002741575.1_ASM274157v1_]  MULTISPECIES: alpha-D-ribose'
                            + ' 1-methylphosphonate 5-phosphate C-P-lyase [Enterobacteriaceae]'
                            + ' [start:4548948] [end:4549793] [genome:GCF_002741575.1_ASM274157v1_]</cite>'))

                tags.li(raw('<i><b>genes_info.tsv</b></i> — tab separated table with 2 columns: <b>gene</b> and'
                            + ' <b>info</b>. Contains information for each gene in the proteome. Information'
                            + ' is extracted from fasta header.'))

                tags.li(raw('<i><b>cluster_info.json</b></i> — contains counted information about the genes in'
                            + ' each cluster. You can use it to take a look at what genes the cluster consists of'
                            + ', and to set a name for the cluster. Here is the example of visualisation'
                            + ' of this information for 1 cluster:'))

                tags.img(src='../images/cluster_names.png')

                tags.li(raw('<i><b>total_list.csv</b></i> — table with clusters in rows and genomes in columns. Number'
                            + ' on the intersection means, how many genes from this cluster the genome contains.'))

                with tags.li(raw('<i><b>clusters.csv</b></i> — data frame with following columns:')):
                    with tags.ol():
                        tags.li(raw('<i>cluster</i> — cluster of a sequence'))
                        tags.li(raw('<i>representative</i> — boolean variable, which indicates, if the given'
                                    + ' sequence is representative'))
                        tags.li(raw('<i>gene</i> — name of the sequence from fasta header with appended genome ID'))
                        tags.li(raw('<i>length</i> — length of the gene\'s sequence'))
                        tags.li(raw('<i>identity</i> — identity of a sequence to representative sequence'
                                    + ' of its cluster'))
                        tags.li(raw('<i>name</i> — name of the sequence from fasta header without appended genome ID'))
                        tags.li(raw('<i>genome</i> — genome id from fasta header'))
                        tags.li(raw('<i>info</i> — information about gene from its fasta header'))

            tags.p(raw('Visualisation of the total list of genes. Each row represents a genome. Each column'
                       + ' represents a cluster (only clusters with <b>size > 5</b> included).'
                       + ' If the genome has the gene, there is a white dot on the intersection.'
                       + ' Otherwise — a black dot.'))

            tags.img(src='../plots/genes_presence.png', alt='genes presence', cls='big_img')

            tags.p(raw('Clustering of genes. These illustrations are for you to take a look at the data.'
                       + ' For higher quality illustrations, and the following analysis it is recommended to perform '
                       + 'these operations in <a href="https://software.broadinstitute.org/morpheus/">Morpheus tool</a>'
                       + ' (also avialable as <a href="https://github.com/cmap/morpheus.R">R package</a>).'))

            tags.img(src='../plots/genes_presence_clustering.png', alt='genes presence clustering')

            tags.p(raw('If you had any problems during the work of the program, or have any notes, please contact'
                       + ' <a href="mailto:vladimirshitov98@gmail.com">vladimirshitov98@gmail.com</a>'))

    with open('./results/summary.html', 'w') as f:
        f.write(str(doc))
