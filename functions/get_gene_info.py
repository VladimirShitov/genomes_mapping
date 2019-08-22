from Bio import Entrez

import constants

Entrez.email = constants.EMAIL


def get_gene_info(gene_id):
    """
    Get gene short name and information by given id.

    Perform NCBI gene search twice: for a given id (it is a locus tag, probably) and for UID found. If given
    gene_id is not found in 'gene' database, performs search in 'protein' database.

    Parameters
    ----------
    gene_id : string to search in NCBI gene database.

    Returns
    -------
    name : str
        Attribute 'Name' from summary for found gene. Usually a short name like 'ytjA'. None if not found.
    info : str
        Attribute 'Summary' from summary for found gene. Brief description of a gene. 'Not found' if not found.

    Examples
    --------
    >>> name, info = get_gene_info('NP_414542.1')
    >>> print(name)
    thrL
    >>> print(info)
    The ThrL leader peptide controls by attenuation the expression of the  thrLABC operon,
    which encodes four out of the five enzymes of threonine biosynthesis pathway,
    in response to the threonine and isoleucine levels . [More information is available at EcoCyc: EG11277].
    """

    db = 'gene'

    search = Entrez.esearch(db, 'Escherichia coli[orgn] {}'.format(gene_id))
    record = Entrez.read(search)

    if not record['IdList']:
        db = 'protein'

        search = Entrez.esearch(db, 'Escherichia coli[orgn] {}'.format(gene_id))
        record = Entrez.read(search)

        if not record['IdList']:
            return None, 'Not found'

    gene_id = record['IdList'][0]

    handle = Entrez.esummary(db=db, id=gene_id)
    record = Entrez.read(handle)

    if db == 'gene':
        summary_dict = record['DocumentSummarySet']['DocumentSummary'][0]
        name = summary_dict['Name']
        info = summary_dict['Summary']

    elif db == 'protein':
        summary_dict = record[0]
        name = summary_dict['Title']
        info = summary_dict['Comment']

    return name, info
