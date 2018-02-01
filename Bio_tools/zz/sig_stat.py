import pandas as pd

for sample in ['stim_WT_CD8', 'stim_DKO_CD8']:
	sig_df = pd.read_csv('%s.spline_pass2.significances.txt.gz' % sample, header=0, sep='\t')
	print sample
	print '0.05', sig_df[sig_df['q-value'] < 0.05].shape[0]
	print '1e-4', sig_df[sig_df['q-value'] < 1e-4].shape[0]
	print '1e-10', sig_df[sig_df['q-value'] < 1e-10].shape[0]
	print
