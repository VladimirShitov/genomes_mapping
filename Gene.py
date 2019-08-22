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
        Nucleotide or aminoacid sequence of a gene.

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

    def __repr__(self):
        """Return `self.get_fasta()` output."""
        return self.get_fasta()

    def get_fasta(self):
        """Return string with a gene in fasta format.

        Returns
        -------
        fasta : str
            String with a gene in fasta format. Format is modified so that Gene() object could be completely
            restored from it. In general format looks like this:
            ">`self.id` `self.info` [start:`self.start`] [end:`self.end`] [genome: `self.genome`]
            `self.sequence`"
        """
        return ('>'
                + self.id
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
