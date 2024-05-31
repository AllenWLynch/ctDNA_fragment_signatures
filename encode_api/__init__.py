from .encode import FeatureConfig


def generate_jinja_config(template, *args, **kwargs):
        
    from jinja2 import Environment, FileSystemLoader
    import os

    environment = Environment(loader=FileSystemLoader(
        os.path.join(os.path.dirname(__file__), "templates/")
    ))
    template = environment.get_template(template)
    
    return template.render(*args, **kwargs)



def make_fragments_config(tissue_type, output, 
                          cell_line=False, exact_ontology=False):
        
    collect=["transcription","H3K4me3", "H3K4me1", "H3K27ac", 
            "H3K9ac", "H3K36me3", "H3K27me3",
            "H3K9me3", "DNase", "POLR2A", "CTCF"]

    features = FeatureConfig(
            tissue_type,
            *collect,
            exact_ontology=exact_ontology,
            cell_line=cell_line,
        ).parse_url(
            genome_build='GRCh38',
        )
    
    template=generate_jinja_config(
        'fragments-template.yaml',
        features,
    )

    with open(output, 'w') as f:
        f.write(template)
