#!/usr/share/gemini/anaconda/bin/python -E
# -*- coding: utf-8 -*-
# usage: gemini_summarize.py <query> <gemini.db>

import sys
import locale
from gemini import GeminiQuery

DP_THRESHOLD = 8
GQ_THRESHOLD = 20

reload(sys)
sys.setdefaultencoding(locale.getpreferredencoding())

gq = GeminiQuery(sys.argv[2], include_gt_cols=True)
gq.run(sys.argv[1], None, True)

header_printed = False
genotype_columns = [ 'gt_depths', 'gt_ref_depths', 'gt_alt_depths', 'gts', 'gt_quals' ]
for row in gq:
	columns = row.row.keys()[:-11]
	if (not header_printed):
		# Print file header
		print '\t'.join(['SAMPLE_ID'] + [ s[:-1] for s in genotype_columns ] + columns)
		header_printed = True

	# Output only het & hom alt variants
	for sample in row['variant_samples']:
		# Drop low depth and low quality variants
		if row['gt_depths'][gq.sample2index[sample]] < DP_THRESHOLD or row['gt_quals'][gq.sample2index[sample]] < GQ_THRESHOLD: continue

		# Print selected genotype and varinat columns
		genotype = [ str(row[col][gq.sample2index[sample]]) for col in genotype_columns ]
		print '\t'.join([sample] + genotype + [ str(row[col]) for col in columns ])
