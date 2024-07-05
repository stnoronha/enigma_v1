# for accessing enrichr API after significant genes are obtained

import json
import requests
from gene_expression import store_list

genes_lst = [
'GATB',
'RAB3C',
'KLF6',
'CYP51A1',
'HDAC9',
'ADCY7',
'FREM3',
'ACTC1',
'FBXL2',
'PCDH10',
'VIT',
'ADTRP',
'SPON2',
'RTP1'
]

def enrichr(gene_lst:list, description = 'Example gene list'):
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(gene_lst)
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)

    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = data['userListId']
    gene_set_library = 'KEGG_2015'
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
    )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text)
    store_list(data,description)
    

