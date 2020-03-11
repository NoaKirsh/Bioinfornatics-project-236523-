from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import pandas as pd


#Coverage vs Replicates.
def samples_vs_coverage_2d(master_table, title): #2d
    fig = plt.figure()
    plt.scatter(master_table['Num of samples'],
               master_table['Coverage percentage'],
               c=master_table['Score'],
               cmap=cm.gist_rainbow_r,
               s=20)
    plt.title(title)
    plt.xlabel('Num of samples')
    plt.ylabel('Coverage percentage')
    plt.legend()

    colmap = cm.ScalarMappable(cmap=cm.gist_rainbow_r)
    colmap.set_array(master_table['Score'])
    plt.clim(0.5, 1)
    colmap.set_clim(vmin=0.5, vmax=1)
    cb = fig.colorbar(colmap)
    cb.ax.set_title('Score')
    plt.show()


# number of samples vs number of individuals.
def mice_sample_3d(master_table):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(180, 270)
    ax.scatter(master_table['Num of samples'],
               master_table['Samples-mice ratio'],
               master_table['Num of mice'],
               c=master_table['Score'],
               cmap=cm.gist_rainbow_r, s=20)
    ax.set_xlabel('Num of samples')
    ax.set_ylabel('')
    plt.yticks([],[])
    ax.set_zlabel('Num of mice')
    colmap = cm.ScalarMappable(cmap=cm.gist_rainbow_r)
    colmap.set_array(master_table['Score'])
    colmap.set_clim(vmin=0.55, vmax=1)
    cb = fig.colorbar(colmap)
    cb.ax.set_title('Score')
    plt.savefig("master graph.png")
    plt.show()


# 3d  number of samples vs number of individuals vs number of reads.
def micesample(master_table): #3d

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(50, 310)
    ax.scatter(master_table['Num of mice'],
               master_table['Num of samples'],
               master_table['Num of reads'],
               c=master_table['Score'],
               cmap=cm.gist_rainbow_r, s=20)
    # plt.xscale('linear')
    plt.ylim(0, 500)
    ax.set_xlabel('Num of mice')
    ax.set_ylabel('Num of samples')
    ax.set_zlabel('Num of reads')
    # fig.colorbar(ax)
    colmap = cm.ScalarMappable(cmap=cm.gist_rainbow_r)
    colmap.set_array(master_table['Score'])
    cb = fig.colorbar(colmap)
    cb.ax.set_ylabel('Score')
    plt.savefig("mice sample graph.png")
    # plt.title("mice-sample ratio")
    plt.show()


#  number of samples vs number of individuals
def samples_vs_mice(master_table): #2d
    fig = plt.figure()
    plt.scatter(master_table['Num of mice'],
               master_table['Num of samples'],
               c=master_table['Score'],
               cmap=cm.gist_rainbow_r,
               s=20)
    plt.title("Number of individuals VS number of samples")
    plt.xlabel('Num of mice')
    plt.ylabel('Num of samples')
    plt.legend()

    colmap = cm.ScalarMappable(cmap=cm.gist_rainbow_r)
    colmap.set_array(master_table['Score'])
    cb = fig.colorbar(colmap)
    cb.ax.set_title('Score')
    plt.show()

def main():

    ####################################################################
                ########## Uploading the Table ###########
    ####################################################################

    master_table = pd.read_csv('MasterTable.csv', sep=',', lineterminator='\r')
    master_table['Num of reads'].replace(r'\s+|\\n', ' ', regex=True, inplace=True)
    master_table.drop(master_table.tail(1).index, inplace=True)  # drop last n rows

    ####################################################################
                    ########## Graphs ###########
    ####################################################################


    samples_vs_coverage_2d(master_table, "")
    best = master_table[master_table['Score'] > 0.85]
    micesample(master_table)

    budget1 = master_table[master_table['Num of reads'] < 400]
    samples_vs_coverage_2d(budget1, "")

    budget4 = master_table[master_table['Num of reads'] > 400]
    samples_vs_coverage_2d(budget4, "")

    mice_sample_3d(master_table)

if __name__ == '__main__':
    main()