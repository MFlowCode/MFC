#!/usr/bin/env python3

import os

#get .bib for files with arXiv identifier
with open(r'_data/publist.yml') as file:

    for pub in pubs:
        for k,v in pub.items():
            if (k=='url'):
                url=v
            if (k=='arxiv'):
                print('working on ArXiv file ' + url)
                os.system("./scripts/arxiv2bib.py " + str(v) + " > papers/" + url +".txt")

