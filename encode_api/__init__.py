from .encode import FeatureConfig

def generate_jinja_config(template, *args, **kwargs):
        
        from jinja2 import Environment, FileSystemLoader
        import os

        environment = Environment(loader=FileSystemLoader(
                os.path.join(os.path.dirname(__file__), "templates/")
        ))
        template = environment.get_template(template)
        
        return template.render(*args, **kwargs)