import os


class BlastDatabase:

    def __init__(self):
        self.genes = []
        self.written_to_a_file = 0

    def feed_genes(self, genes):
        self.genes.extend(genes)

    def create_db(self, append=False, path='.'):
        flag = 'a' if append else 'w'

        with open(path + 'current_db.faa', flag) as f:
            for gene in self.genes[self.written_to_a_file:]:
                f.write(str(gene))

        self.written_to_a_file += len(self.genes[self.written_to_a_file:])

        os.system('makeblastdb -dbtype prot -in {path}current_db.faa -out {path}cur_db'.format(path=path))
