import csv
import glob
import numpy

import Morph
import Morph.features


data = {}
for annotation in glob.glob('/diskmnt/Users2/avisani/analysis/mCRC/vascular_andre/2D_5K/*s.csv'):
    sample = annotation.split('/')[-1].split('_')[0]
    for cells in glob.glob('/diskmnt/primary/Xenium/data/*/*' + sample + '*/cells.csv.gz'):
        data[sample] = Morph.readers.cells(cells)
        group = {}
        with open(annotation, 'r') as f:
            dict_reader = csv.DictReader(f)
            for line in dict_reader:
                group[line['cell_id']] = line['group']
        barcodes = group.keys()
        types = []
        for g in data[sample]['g']:
            types.append(group[g] if g in barcodes else 'unmatched')
        data[sample]['g'] = types
d = 10
cell_types = set()
for sample in data:
    cell_types = cell_types | set(data[sample]['g'])
cell_types = cell_types - {'T_NK', 'unmatched'}
distance = Morph.features.Distance()
distances = {}
image = {}
total = {}
for sample in data:
    image[sample] = {}
    distances[sample] = {}
    for t in cell_types:
        G = {t}
        image[sample][t] = Morph.backbone(data[sample], ['xenium', d], ['total', G], ['maximum'], ['naive'], ['naive'], ['naive'], ['naive'])
    total[sample] = numpy.zeros_like(image[sample][t])
    for t in cell_types - {'Angiogenic_EC', 'LSEC', 'Lymphatic_EC'}:
        total[sample] += image[sample][t]
    for t in {'Angiogenic_EC', 'LSEC', 'Lymphatic_EC'}:
        distances[sample][t] = distance.maximum(image[sample][t])
cell_types_reference = {'Angiogenic_EC', 'LSEC'}
distribuitions = {}
for t in cell_types - {'Angiogenic_EC', 'LSEC', 'Lymphatic_EC'}:
    distribuitions[t] = {r: {} for r in cell_types_reference}
    for c in cell_types_reference:
        distribuitions[t][c] = {b: {} for b in range(d, 100, d)}
        if t == c:
            continue
        for l, b in zip(range(1, int(100 / d)), range(d, 100, d)):
            for sample in data:
                distribuitions[t][c][b][sample] = sum(image[sample][t][distances[sample][c].astype(int) == l]) / sum(total[sample][distances[sample][c].astype(int) == l]) if sum(total[sample][distances[sample][c].astype(int) == l]) > 0 else 0
with open('adrienne2compare.csv', 'w') as f:
    dict_writer = csv.DictWriter(f, fieldnames=['to', 'from', 'bin'] + [sample for sample in data])
    dict_writer.writeheader()
    for t in cell_types - {'Angiogenic_EC', 'LSEC', 'Lymphatic_EC'}:
        for c in cell_types_reference:
            if t == c:
                continue
            for b in range(d, 100, d):
                row = {sample: distribuitions[t][c][b][sample] for sample in data}
                row['to'] = t
                row['from'] = c
                row['bin'] = b
                dict_writer.writerow(row)
