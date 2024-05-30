import deeptools.getScorePerBigWigBin as score_bw
from deeptools.correlation import Correlation
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import requests
from urllib.parse import urlparse
from collections import defaultdict
import logging

histone_marks = ["H3K4me3", "H3K4me2", "H3K4me1", "H3K27ac", "H3K9ac", "H3K36me3",
                 "H3K27me3", "H3K9me3"]
TF_marks = ['CTCF', 'POLR2A']
DNase = ['DNase']
gene = ['transcription']

logger = logging.getLogger(__name__)
class FeatureConfig:
    def __init__(self, tissue_type, exact_ontology=True, cell_line=False, *featuretypes):
        if exact_ontology:
            self.tissue_type=tissue_type
            self.searchTerm=None
        else:
            self.tissue_type=None
            self.searchTerm=tissue_type

        self.cell_line = cell_line
        self.featureLst = list(featuretypes)          
        
    def find_ENCODE_data_entries(self, assay_type, target=None, tissue_type=None, searchTerm=None, cell_line=False, genome_build="GRCh38", rfa="ENCODE4"):
        if 'histone' in assay_type.lower():
            assay_title = 'Histone+ChIP-seq'
            assay_filter_line = f'&assay_title={assay_title}'
            file_type='bigWig'
            output_type='signal+p-value'

            
    
        elif 'TF' in assay_type.upper():
            assay_title='TF+ChIP-seq'
            assay_filter_line = f'&assay_title={assay_title}'
            if target =="CTCF":
                file_type='bed+narrowPeak'
                output_type='*'
            # RNA polymerase II 
            elif target=="POLR2A":
                file_type='bigWig'
                output_type='signal+p-value'
                

        elif 'RNA' in assay_type.upper():
            assay_titles =[ 'total+RNA-seq', 'polyA+plus+RNA-seq'] # add compatible RNA-seq here
            assay_filter_line = ''
            for assay_title in assay_titles:
                assay_filter_line += f'&assay_title={assay_title}'

            file_type='tsv'
            output_type='gene+quantifications'

        elif 'dnase' in assay_type.lower():
            assay_title='DNase-seq'
            assay_filter_line = f'&assay_title={assay_title}'
            file_type='bigWig'
            output_type='*'
            
    
        if not tissue_type is None:
            tissue_type = tissue_type.replace(' ', '+')
            term_filter_line = f'&biosample_ontology.term_name={tissue_type}' 
        else:
            searchTerm = searchTerm.replace(' ', '+')
            term_filter_line = f'&searchTerm={searchTerm}' 

        
        sample_class = 'tissue' if not cell_line else 'cell+line'

        target_filter_line = f'&target.label={target}' if not target is None else ''
        
        url = (
            'https://www.encodeproject.org/search/?type=Experiment'
            '&status=released'
            f'{assay_filter_line}'
            f'&biosample_ontology.classification={sample_class}'# tissue or cell+line or primary+cell
            '&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens'
            '&replicates.library.biosample.perturbed=false'
            f'{term_filter_line}' 
            '&award.project=ENCODE'
            # '&award.rfa=ENCODE4'
            f'{target_filter_line}'
            '&frame=embedded'
            '&format=json'
        )
        print(url)
        r = requests.get(url)
        experiments = r.json()['@graph']
    
        assert len(experiments)>0,  "No matched biosample entry found!"
            
    
        # Obtain a subset of the biosample metadata.
        biosample_ontologies = []
        for e in experiments:
            # Get the name of the cell line, as well as more
            # general categorical information.
    
            biosample_ontology_condensed = {
                k: v
                for k, v in e.get('biosample_ontology', {}).items()
                if k == 'term_name' or k.endswith('_slims') or k=='summary'
            }
            biosample_ontology_condensed.update(
                {'summary': e.get('simple_biosample_summary')}
            )
            biosample_ontology_condensed.update(
                {'dataset': e.get('@id')}
            )

            biosample_ontology_condensed.update(
                {'assay':e.get('assay_title')}
            )
            biosample_ontologies.append(biosample_ontology_condensed)
        biosample_metadata = pd.DataFrame(biosample_ontologies)
    
        display(biosample_metadata)

        use_rows = list(map(int, 
        input("\nEnter the indices of rows to use : ").strip().split()))
        experiments = [experiments[x] for x in use_rows]
    
        # Parse the experiments into query parameters.
        datasets = [
            '&dataset={}'.format(e.get('@id'))
            for e in experiments
        ]
        url = (
            'https://www.encodeproject.org/search/?type=File'
            '&status=released'
            f'&file_type={file_type}'
            f'&output_type={output_type}'
            f'&assembly={genome_build}'
            f'&award.rfa={rfa}'
            '&field=output_type' 
            '&field=target.label' 
            '&field=award.rfa' 
            '&field=dataset'
            '&field=biological_replicates'
            '&field=cloud_metadata.url'
            '&frame=object'
            '&format=json'
            '&limit=all'
            '&{}'.format(''.join(datasets))
        )
        r = requests.get(url)
    
        # Filter the files to only those belonging to multiple replicates.
        files = [
            f
            for f in r.json()['@graph']
            if len(f['biological_replicates']) >= 1
        ]
    
        # if len(files)==0:
        #     files = [
        #     f
        #     for f in r.json()['@graph']
        #     if len(f['biological_replicates']) == 1
        # ]
        assert len(files)>0, "No matched files found"
    
        file_metadata = pd.DataFrame(files)
        # Flatten JSON.
        try:
            file_metadata['target'] = file_metadata.target.apply(lambda x: x.get('label'))
        except:
            pass
        file_metadata['cloud_metadata'] = file_metadata.cloud_metadata.apply(lambda x: x.get('url'))
        file_metadata['project'] = file_metadata.award.apply(lambda x: x.get('rfa'))
        file_metadata['output_type'] = file_metadata.output_type
    
        # Rename column.
        file_metadata = file_metadata.rename(columns={'@id': 'file'})
        # Join with biosample_metadata.
        merged_df = file_metadata.merge(biosample_metadata, how='outer', on='dataset')
        column_order = [
            'file',
            'dataset',
            'project',
            'biological_replicates',
            'output_type',
            'target',
            'cell_slims',
            'developmental_slims',
            'summary',
            'organ_slims',
            'system_slims',
            'term_name',
            'cloud_metadata',

        ]
        merged_df = merged_df[[col for col in column_order if col in merged_df.columns]]
        merged_df = merged_df.dropna().sort_values(by=['dataset', 'project']).reset_index(drop=True)

        return merged_df
    
    def select_files(self, assay_type, feature, genome_build, rfa):
        try:
            merged_df = self.find_ENCODE_data_entries(assay_type, target=feature, tissue_type=self.tissue_type, searchTerm=self.searchTerm, cell_line=self.cell_line,\
                                                      genome_build=genome_build, rfa=rfa)
        except AssertionError as e :
            logger.warning(e)   
            return None
        
        display(merged_df)    
        use_rows = list(map(int, 
        input("\nEnter the indices of rows to use : ").strip().split()))    
        return merged_df.iloc[use_rows]

    def parse_url(self,  genome_build="GRCh38", rfa="ENCODE4"):
        if genome_build=='hg19':
            rfa='ENCODE3' # if hg19 provided as genome build version, rfa should be overwritten to ENCDOE3
        feature_url_dict=defaultdict(list)
        for feature in self.featureLst:
            if feature in histone_marks:
                assay_type="histone chip-seq"
                final_df = self.select_files(assay_type, feature, genome_build, rfa)
                
                if final_df is None:
                    continue
                for url in final_df.cloud_metadata.values:
                    feature_url_dict[feature].append(url)


            elif feature in TF_marks:
                assay_type="TF chip-seq"
                final_df = self.select_files(assay_type, feature, genome_build, rfa)
                
                if final_df is None:
                    continue
                for url in final_df.cloud_metadata.values:
                    feature_url_dict[feature].append(url)

            elif feature in gene:
                assay_type="RNA-seq"
                final_df = self.select_files(assay_type, None, genome_build, rfa)
                
                if final_df is None:
                    continue
                for url in final_df.cloud_metadata.values:
                    feature_url_dict[feature].append(url)     

            elif feature in DNase:
                assay_type="DNase-seq"
                final_df = self.select_files(assay_type, None, genome_build, rfa)
                if final_df is None:
                    continue
                for url in final_df.cloud_metadata.values:
                    feature_url_dict[feature].append(url)  
                
                
                    
            else:
                logger.info(f"{feature} is currently not supported!")
                

        return feature_url_dict


    def generate_draft_config(self):
        pass

    
