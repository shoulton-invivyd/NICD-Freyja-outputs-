import os


with open('recent_samples.txt') as file:
    samples = [line.rstrip() for line in file]


vars_orig = [var for var in os.listdir('../variants/')]
vars = [var.split('_S')[0].split('.tsv')[0].replace('_','-') for var in os.listdir('../variants/')]
    

recent_inds = [j for j,v in enumerate(vars)  if v in samples]

recent_samples = ['../variants/'+vars_orig[j] for j in recent_inds]

for r in recent_samples:
    print(r)

# with open('recent_samples.txt') as file:
#     samples = [line.rstrip() for line in file]