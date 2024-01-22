import sys, os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from os import path, walk
from sklearn.decomposition import PCA
import itertools

from Afanc.utilities.runCommands import command


def makePCA(distance_matrix, species_dict):
    # Perform PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(distance_matrix.values)

    # Create a DataFrame for visualization
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
    pca_df.index = distance_matrix.index  # Assuming index contains sequence IDs

    # Map species to the DataFrame based on the dictionary
    pca_df['Species'] = pca_df.index.map(species_dict)

    # Create a classifier dictionary
    unique_species = pca_df['Species'].unique()
    # classifier_dict = dict(zip(unique_species, plt.cm.tab10(np.arange(len(unique_species)))))
    markers = list(itertools.product(['o', 's', '^', 'D', 'v'], plt.cm.tab10.colors))
    classifier_dict = dict(zip(unique_species, markers))

    # Plot the PCA with colored points
    fig, ax = plt.subplots()
    for species, (shape, color) in classifier_dict.items():
        subset_df = pca_df[pca_df['Species'] == species]

        # # Calculate the mode for each cluster
        # for cluster in subset_df.groupby('Species').groups.values():
        #     # Filter the DataFrame using boolean indexing
        #     cluster_data = subset_df.loc[cluster]
            
        #     # Calculate the mode
        #     cluster_mode = cluster_data[['PC1', 'PC2']].mode().iloc[0]
            
        #     # Draw a ring around the mode with a radius of 0.1 times the range
        #     radius = 0.1 * np.ptp(cluster_data[['PC1', 'PC2']].values)
        #     circle = plt.Circle((cluster_mode['PC1'], cluster_mode['PC2']), radius, color='black', fill=False)
        #     ax.add_patch(circle)

        ax.scatter(subset_df['PC1'], subset_df['PC2'], marker=shape, color=color, label=species)

    ax.set_xlabel('Principal Component 1')
    ax.set_ylabel('Principal Component 2')
    ax.grid(True)
    ax.legend()

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ## add a scree plot
    # scree_ax = fig.add_axes([0.69, 0.67, 0.2, 0.2])
    # explained_variance_ratio = pca.explained_variance_ratio_
    # scree_ax.plot(range(1, len(explained_variance_ratio) + 1), np.cumsum(explained_variance_ratio), marker='o', linestyle='--')
    # scree_ax.set_xlabel('N. PCs')
    # scree_ax.set_ylabel('Cum. Variance')

    # scree_ax.axvline(x=2, color='red', linestyle='--')  ## Draw a vertical line at x=2

    # y_intersect = np.cumsum(explained_variance_ratio)[1]
    # scree_ax.text(2, y_intersect, f'(2,{y_intersect:.2f})', color='red', verticalalignment='bottom')

    # scree_ax.grid(True)

    plt.show()



def readDistOut(dist_file, id_dict):

    # Read the file into a DataFrame
    df = pd.read_csv(dist_file, sep='\t', header=None, names=["ref_path", "query_path", "mash_dist", "p", "matching_hashes"])

    # Extract file names from paths
    df['ref_ID'] = df['ref_path'].apply(lambda x: path.basename(x))
    df['query_ID'] = df['query_path'].apply(lambda x: path.basename(x))

    # Create a pivot table to construct the distance matrix
    distance_matrix = df.pivot(index='ref_ID', columns='query_ID', values='mash_dist')

    # Fill the diagonal with zeros
    distance_matrix = distance_matrix.fillna(0)

    # Fill in the missing values by mirroring the existing values
    distance_matrix = distance_matrix + distance_matrix.T

    # Save the result to a CSV file
    distance_matrix.to_csv("mash_out.dist")

    return distance_matrix


def flattenDirectory(directory_path):
    flat_dir_dict = {}

    for root, dirs, files in os.walk(directory_path):
        # Extract the subdirectory name from the current path
        subdirectory = path.relpath(root, directory_path)

        # Add non-directory files to the result dictionary
        files_in_subdirectory = [path.join(root, f) for f in files if os.path.isfile(os.path.join(root, f))]
        flat_dir_dict[subdirectory] = files_in_subdirectory

    return flat_dir_dict


def mash(prefix, fastas):
    """ run mash
    """
    
    mashdist_out = path.abspath(f"{prefix}_mashdist.txt")

    mash_sketchline = f"mash sketch -o ref {' '.join(fastas)}"
    mash_distline = f"mash dist ref.msh ref.msh > {mashdist_out}"
    command(mash_sketchline, "MASH").run_comm(0, None, None)
    command(mash_distline, "MASH").run_comm(0, None, None)

    return mashdist_out

def fastANI(prefix, fastas):
    """ run fastANI
    """

    outfile = f"{prefix}_fastANI.txt"

    with open(f"fastas.txt", 'w') as fout:
        print('\n'.join(fastas), file=fout)

    argline = f"fastANI --ql fastas.txt --rl fastas.txt -o {outfile} -t 10"
    print(argline)
    # subprocess.call(argline, shell=True)

    return outfile

if __name__=="__main__":
    fasta_dir="/home/amorris/BioInf/afanc_kleb_WD/klebsiella_3.0_getds/"
    out_dir="/home/amorris/BioInf/afanc_kleb_WD/"

    taxon_dict = flattenDirectory(fasta_dir)
    reversed_dict = {path.basename(value): key for key, values in taxon_dict.items() for value in values}

    fasta_list = [fasta for fastas in taxon_dict.values() for fasta in fastas]
    dist_out = mash("kleb", fasta_list)
    # dist_out = fastANI("kleb", fasta_list)

    distance_matrix = readDistOut(dist_out, reversed_dict)

    makePCA(distance_matrix, reversed_dict)