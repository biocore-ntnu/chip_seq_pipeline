
# from jinja2 import Environment, FileSystemLoader

# env = Environment(loader=FileSystemLoader("utils/"))
# env.get_template("config_template.yaml").render()

# have list of all config files
# read their values
# if the values do not exist, replace with default

import yaml

template_path = "utils/config_description.yaml"

config_files = ["config.yaml", "tests/test_data/paired_end/config.yaml",
                "tests/test_data/dna_repair/config.yaml",
                "tests/test_data/keep_the_tips/config.yaml",
                "example/config.yaml"]

template = yaml.load(template_path)


for config_file in config_files:

    # Warn about values not there
    # Warn about values there, but not in template
    # Create backup?

    config = yaml.load(open(config_file))

    # for k, v in template():
    #     if k in template
