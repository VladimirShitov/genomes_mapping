class Gene:
    """Gene representation in form of a python class.

    Attributes
    ----------
    id : str
        Any name to identificate gene later. Preferably, locus tag.
    info : str
        Any additional information about the gene. Usually part of a fasta file after a gene name, but
        before sequence.
    sequence : str
        Aminoacid sequence of a gene.
    nucleotide_sequence : str
        Nucleotide sequence of a gene.

    Methods
    -------
    get_fasta()
        Return string with a gene in fasta format.
    to_file(path)
        Write gene in a fasta format to the file.
    __repr__
        Return `self.get_fasta()` output.
    """

    def __init__(self, gene_id, info, sequence, genome):
        """
        Parameters
        ----------
        gene_id : str
            Any name to identificate gene later. Usually, part of a fasta header right after '>'.
        info : str
            Any additional information about the gene. Usually part of a fasta file after a gene name, but
            before sequence.
        sequence : str
            Nucleotide or aminoacid sequence of a gene.
        """
        self.id = gene_id
        self.info = info
        self.sequence = sequence
        self.genome = genome
        self.start = None
        self.end = None
        self.nucleotide_sequence = None

    def __repr__(self):
        """Return `self.get_fasta()` output."""
        return self.get_fasta()

    def get_fasta(self, join_genome_to_name=False):
        """Return string with a gene in fasta format.

        Parameters
        ----------
        join_genome_to_name : bool
            If True, append [genome: `self.genome`] to gene id in fasta header. It may be useful if you need to
            differ 2 genes with the same name but from different genomes.

        Returns
        -------
        fasta : str
            String with a gene in fasta format. Format is modified so that Gene() object could be completely
            restored from it. In general format looks like this:
            ">`self.id` `self.info` [start:`self.start`] [end:`self.end`] [genome: `self.genome`]
            `self.sequence`"
        """
        gene_id = self.id
        if join_genome_to_name:
            gene_id += '[genome:{}] '.format(self.genome)

        return ('>'
                + gene_id
                + ' '
                + self.info
                + ' '
                + '[start:{}] '.format(self.start)
                + '[end:{}] '.format(self.end)
                + '[genome:{}] '.format(self.genome)
                + '\n'
                + self.sequence
                + '\n'
                )

    def to_file(self, path):
        """Write gene in a fasta format to the file.

        Parameters
        ----------
        path : str
            Path to a file, where to write gene.
        """
        with open(path + self.id, 'w') as f:
            f.write(self.get_fasta())
