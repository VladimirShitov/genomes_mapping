import pandas as pd
from Gene import Gene


class Genome:
    """Genome files, information and methods for analysis.

    Attributes
    ----------
    files : dict
        Dictionary with keys containing names of all possible files avialable for a genome on NCBI FTP-site.
        https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/ chapter 11. Values contain path to this file
        on a local computer or None.
    raised_errors : bool
        False by default, True if any error occured while working with a genome.
    log : str
        Log with all information about errors.
    folder_path : str
        Path to a folder, where files for genome can be found.
    files_prefix : str
        Prefix for a files in folder. Usually contains assembly number.
        I.e. for a path '/home/user/E_coli/GCF_000007445.1_ASM744v1_protein.faa' `folder_path` would be
        '/home/user/E_coli/' and `files_prefix` - 'GCF_000007445.1_ASM744v1_'
        TODO: delete this attribute. Can check file_name.endswith('protein.faa')
    feature_table : pandas.DataFrame
        feature_table.txt file in dataframe format from pandas after self.read_feature_table() method called.
    genes : list
        List of genes from protein.faa file for a genome. Each gene is an instance of Gene() class.
        Filled after self.read_protein_faa() methid called.

    Methods
    -------
    read_feature_table()
        Read feature_table.txt in a `folder_path`. Save pandas.DataFrame in `self.feature_table`.
    read_protein_faa()
        Read protein.faa in a `folder_path`. Save set of Gene() objects in `self.genes`.
    __repr__
        Return `folder_path` and `files_prefix` of a genome
    get_gene_by_id(gene_id)
         Get gene, which `id` is the same as `gene_id` parameter.
    set_gene_positions()
         Set `start` and `end` for each gene in `self.genes`. Requires `self.feature_table` read.
    write_genes_to_file(path)
        Create a file and write all genes sequences in fasta format in it.

    """

    def __init__(self, folder_path, files_prefix):
        """Set genome folder path and prefix for files.

        Parameters
        ----------
        folder_path : str
            Path to a folder, where files for genome can be found.
        files_prefix : str
            Prefix for a files in folder. Usually contains assembly number.
            I.e. for a path '~/Documents/E_coli/GCF_000007445.1_ASM744v1_protein.faa' `folder_path` would be
            '~/Documents/E_coli/' and `files_prefix` - 'GCF_000007445.1_ASM744v1_'
        """
        self.files = {'assembly_report.txt': None,
                      'assembly_stats.txt': None,
                      'assembly_regions.txt': None,
                      'assembly_structure directory': None,
                      'cds_from_genomic.fna': None,
                      'feature_count.txt': None,
                      'feature_table.txt': None,
                      'genomic.fna': None,
                      'genomic.gbff': None,
                      'genomic.gff': None,
                      'genomic.gtf': None,
                      'genomic_gaps.txt': None,
                      'protein.faa': None,
                      'protein.gpff': None,
                      'rm.out': None,
                      'rm.run': None,
                      'rna.fna': None,
                      'rna.gbff': None,
                      'rna_from_genomic.fna': None,
                      'translated_cds.faa': None,
                      'wgsmaster.gbff': None,
                      }

        self.raised_errors = False
        self.log = ''

        self.folder_path = folder_path
        self.files_prefix = files_prefix
        self.genes = []
        self.feature_table = None

        for key in self.files:
            self.files[key] = self.folder_path + self.files_prefix + key

    def __repr__(self):
        """Return `folder_path` and `files_prefix` of a genome"""
        return 'Path: {} Prefix: {}'.format(self.folder_path, self.files_prefix)

    def __getitem__(self, index):
        """Get gene by it's name through an index."""
        return self.get_gene_by_id(index)

    def read_feature_table(self, only_protein_coding=True):
        """Read feature_table.txt in a `folder_path`. Save pandas.DataFrame in `self.feature_table`.

        Parameters
        ----------
        only_protein_coding : bool, optional
            If True, feature_table will be filtered as following:
            '# feature' in ['CDS', 'gene', 'operon']
            'class' in ['protein_coding', 'with_protein']
        """
        try:
            feature_table = pd.read_csv(self.files['feature_table.txt'], sep='\t')

            if only_protein_coding:
                feature_table = feature_table[(feature_table['# feature'].isin(['CDS', 'gene', 'operon'])) &
                                              (feature_table['class'].isin(['protein_coding', 'with_protein']))]

            self.feature_table = feature_table.dropna(subset=['product_accession'])

        except FileNotFoundError:
            self.raised_errors = True
            self.log += 'Cannot get feature_table\n'

    def read_protein_faa(self):
        """Read protein.faa in a `self.folder_path`. Save list of Gene() objects in `self.genes`.

        Examples
        --------
        >>> K12_PATH = '~/Documents/E_coli/K12/'
        >>> K12_PREFIX = 'GCF_000005845.2_ASM584v2_'
        >>> k12 = Genome(K12_PATH, K12_PREFIX)
        >>> print(len(k12.genes))
        0
        >>> k12.read_protein_faa()
        >>> print(len(k12.genes))
        4242
        >>> print(k12.genes[0])
        >NP_414542.1 thr operon leader peptide [Escherichia coli str. K-12 substr. MG1655]
        MKRISTTITTTITITTGNGAG

        """
        try:
            with open(self.files['protein.faa'], 'r') as f:
                genes = f.read()
            genes = [x for x in genes.split('>') if x != '']

            for gene in genes:
                info = gene.split('\n')[0]
                gene_id = info.split(' ')[0]
                info = ' '.join(info.split(' ')[1:])

                seq = ''.join(gene.split('\n')[1:-1])  # The last element is always ''

                self.genes.append(Gene(gene_id=gene_id,
                                       info=info,
                                       sequence=seq,
                                       genome=self.files_prefix))

        except FileNotFoundError:
            self.raised_errors = True
            self.log += 'Cannot get protein.faa\n'

    def get_gene_positions(self, gene=None, gene_id=None):
        """Return start and end for a given gene.

        Parameters
        ----------
        gene : Gene
            Instance of an object Gene() for which to return positions. If None, `gene_id` must be given.
            Ignored, if `gene_id` is given.
        gene_id : str
            Particular object of class Gene() `id` attribute

        Returns
        -------
        start
            Integer number with start position of a gene or None if gene or position was not found.
        end
            Integer number with end position of a gene or None if gene or position was not found.

        Examples
        --------
        >>> K12_PATH = '~/Documents/E_coli/K12/'
        >>> K12_PREFIX = 'GCF_000005845.2_ASM584v2_'
        >>> k12 = Genome(K12_PATH, K12_PREFIX)
        >>> k12.read_protein_faa()
        >>> k12.read_feature_table()
        >>> print( k12.get_gene_positions(gene=k12.genes[0]) )
        (190, 255)
        >>> print( k12.get_gene_positions(gene_id='NP_414542.1') )
        (190, 255)
        """
        if gene_id:
            gene = self.get_gene_by_id(gene_id)

        row = self.feature_table[(self.feature_table['product_accession'] == gene.id)]
        row = row.dropna(subset=['start', 'end'])

        if row.shape[0] > 0:
            row = row.iloc[0]
        else:
            return None, None

        start = row['start']
        end = row['end']

        return start, end

    def set_gene_positions(self):
        """Set `start` and `end` for each gene in `self.genes`. Required `self.feature_table` read."""
        for gene in self.genes:
            gene.start, gene.end = self.get_gene_positions(gene=gene)

    def get_gene_by_id(self, gene_id):
        """Get gene, which `id` is the same as `gene_id` parameter."""
        for gene in self.genes:
            if gene.id == gene_id:
                return gene

    def write_genes_to_file(self, path):
        """Create a file and write all genes sequences in fasta format in it."""
        with open(path, 'w') as f:
            for gene in self.genes:
                f.write(gene.get_fasta())
