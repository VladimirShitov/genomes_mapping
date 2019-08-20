def plot_alignment_distribution(alignment_list, only_multiple=False, log_scale=False, **kwargs):
    lengths = np.array([len(alignment_list[x]) for x in alignment_list])
    if only_multiple:
        lengths = lengths[lengths > 1]

    plt.figure(figsize=(9, 5))
    sns.distplot(lengths, kde=False, **kwargs)
    if log_scale:
        plt.yscale('log')

    plt.title('Multiple aligned genes distribution')
    plt.xlabel('Number of alignments to 1 gene')
    plt.ylabel('Number of genes')

# TODO: document
