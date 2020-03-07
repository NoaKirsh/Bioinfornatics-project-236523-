from datetime import datetime
import numpy as np
import pandas as pd
import random
import os

# Chooses random samples for a given number of mice and percentage of samples.
# Attempts to keep the number of Active and Unstimulated cells similar.

def choosing_samples(df, genes, num_of_mice=2, percentage_of_samples=20):
    names_of_mice = df['Individuals'].unique().tolist()  # get names of all individuals
    random_mice = random.sample(names_of_mice, num_of_mice)  # choose randomly the individuals (by name)
    df = df[df['Individuals'].isin(random_mice)]  # filter the df - get only the chosen mice

    num_of_active_to_keep = max(int(len(df[df['Stimulus'] == 'Active']) * (percentage_of_samples/100)), 1)
    num_of_unstimulated_to_keep = max(int(len(df[df['Stimulus'] == 'Unstimulated']) * (percentage_of_samples/100)), 1)
    active_df = df[df['Stimulus'] == 'Active'].copy()
    unstimulated_df = df[df['Stimulus'] == 'Unstimulated'].copy()
    if len(active_df) == 0 or len(unstimulated_df) == 0:
        raise Exception('oppssssss')
    active_to_keep = active_df.sample(n=num_of_active_to_keep, axis=0)
    unstimulated_to_keep = unstimulated_df.sample(n=num_of_unstimulated_to_keep, axis=0)
    df = pd.concat([active_to_keep, unstimulated_to_keep]).sort_values('id')
    temp_genes = genes.drop(axis=1, columns=df['CellName'])
    genes = genes.drop(axis=1, columns=temp_genes.columns.values)
    return genes

# Generates random gene samplings from the data, based on their probability to be chosen in the “real world” for different
# levels of percentage of the data.

def coverage(df_src, percentage=100):
    df = df_src.copy()
    df_prob = df.copy()
    for col in df.columns:
        df[col].values[:] = 0
    sum_samples = df_prob.sum(axis=0).tolist()  # vector of sums of all columns - number of total reads from each sample.
    prob_samples = [x / sum(sum_samples) for x in sum_samples]
    df_prob = df_prob.astype(float)
    for j, col in zip(range(df_prob.shape[1]), df_prob):  # columns
        for i in range(df_prob.shape[0]):
            df_prob.iat[i, j] = df_prob.iloc[i][j] / sum_samples[j]
    col_amounts = np.random.multinomial(n=sum(sum_samples)*(percentage/100), pvals=prob_samples) # how many reads for each sample
    for col in range(0, df.shape[1]):
        gene_amounts = np.random.multinomial(n=col_amounts[col], pvals=df_prob.iloc[:, col]) # how many reads for each gene in the sample
        for row in range(df.shape[0]):
            df.iat[row, col] = gene_amounts[row]
    return df

# Creates csv names with the appropriate data.

def creates_csv_name(csv_type, num_of_good_samples, num_of_mice, percentage_of_samples, coverage_percentage):
    name = str(num_of_good_samples)\
           + '_' + str(num_of_mice) \
           + '_' + str(percentage_of_samples) \
           + '_' + str(coverage_percentage) \
           + '_' + csv_type \
           + '.csv'
    return name


def main():
    ###################-----load data-----########################

    rawcounts = pd.read_csv('raw_data.tab', sep='\t', lineterminator='\r')
    rawcounts['gene'].replace(r'\s+|\\n', ' ', regex=True, inplace=True)
    metadata = pd.read_csv('metadata.tab', sep='\t', lineterminator='\r')
    metadata['id'].replace(r'\s+|\\n', ' ', regex=True, inplace=True)
    ERCC_conc = pd.read_csv('ERCC_conc.tab', sep='\t', lineterminator='\r')
    ERCC_conc['ERCC ID'].replace(r'\s+|\\n', ' ', regex=True, inplace=True)
    ERCC_conc = ERCC_conc.sort_values('ERCC ID')
	
    ###################-----filter ERCC-----########################
	
    metadata.loc[(metadata['Stimulus'] == 'Activated'),'Stimulus'] = 'Active'
    metadata = metadata[metadata['Age'] == 'Young']
    cols_to_keep = metadata['CellName'].tolist()
    cols_to_keep.append('gene')
    temp_rawcounts = rawcounts.drop(axis=1, columns=cols_to_keep)
    rawcounts = rawcounts.drop(axis=1, columns=temp_rawcounts.columns.values)
    del temp_rawcounts
    genes = rawcounts.head(-92)
    genes = genes[genes.sum(axis=1) > 5000.0]
    ERCC = rawcounts.tail(92)
	
    ###################-----ERCC Corrolation-----########################
	# removes all the low quality data based on the correlation to the ERCC genes.
	
    ERCC2 = ERCC.copy()
    del ERCC2['gene']
    corrs = ERCC2.head(1).copy().rename(index={31089: "correlation"}).transpose()
    ERCC3 = ERCC2.assign(concentration=ERCC_conc['concentration'].tolist())
    for sample in ERCC2.columns:
        corrs.loc[sample, 'correlation'] = ERCC3[sample].corr(ERCC3['concentration'])
    del ERCC2, ERCC3
    corrs = corrs.sort_values('correlation', ascending=False)
	
    ###################-----Choosing samples by: *ERCC number *num_of_mice *percentage_of_samples-----########################
	
    for num_of_good_samples in range(300, 533, 20): #every loops different number of best samples available
		best_metadata = metadata[metadata.CellName.isin(corrs.head(num_of_good_samples).index.values)]
		for  in range(1, 10):# 1 to 9
			for percentage_of_samples in [60, 70, 80, 90, 100]:
				chosen_sapmles_data = choosing_samples(best_metadata, genes, num_of_mice, percentage_of_samples)
				for coverage_percentage in range(40, 141, 10):
					print('num_of_good_samples = ' + str(num_of_good_samples) + ', num_of_mice = '+str(num_of_mice) +
						  ', percentage_of_samples = '+str(percentage_of_samples) +
						  ', coverage_percentage '+str(coverage_percentage) + ' ' + str(datetime.now()))
					try:
						chosen_data = coverage(chosen_sapmles_data, coverage_percentage)
					except:
						print('errooooooorrrrr')
						continue
					chosen_data.insert(loc=0, column='gene', value=genes['gene'])
					chosen_metadata = best_metadata[best_metadata.CellName.isin(chosen_data.columns.values)]
					chosen_data_name = creates_csv_name('rawcounts', num_of_good_samples, num_of_mice, percentage_of_samples, coverage_percentage)
					chosen_metadata_name = creates_csv_name('metadata', num_of_good_samples, num_of_mice, percentage_of_samples, coverage_percentage)
					csv_path_rowcounts = os.path.abspath(chosen_data_name)
					csv_path_metadata = os.path.abspath(chosen_metadata_name)
					chosen_data.to_csv(csv_path_rowcounts, index=None, header=True)
					chosen_metadata.to_csv(csv_path_metadata, index=None, header=True)
    print('done')


if __name__ == '__main__':
    print('start '+str(datetime.now()))
    main()
    print('end ' + str(datetime.now()))
