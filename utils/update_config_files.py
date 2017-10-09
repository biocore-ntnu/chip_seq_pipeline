from jinja2 import Environment, FileSystemLoader
import os
import subprocess

from time import ctime
import yaml
from collections import OrderedDict

from textwrap import wrap

def ordered_load(stream, Loader=yaml.Loader, object_pairs_hook=OrderedDict):
    class OrderedLoader(Loader):
        pass
    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))
    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)
    return yaml.load(stream, OrderedLoader)

j2_env = Environment(loader=FileSystemLoader("."),
                    trim_blocks=True)

def format_descriptions(desc):
    capitalized = " ".join([s.capitalize() for s in desc.split("_")])
    dashes = "-" * len(desc)
    return " ".join(["#", dashes, "\n#", capitalized, "\n#", dashes])

def format_comments(c):

    return "# " + "\n# ".join(wrap(str(c), 78))

def required(r):

    return "# (Required)" if r else "# (Not required)"

def default(d):

    if d != "None" and d:
        return str(d)
    else:
        return ""

def example(e):

    if not e:
        return ""

    if isinstance(e, list):
        example = "\n#".join([str(s) for s in e])
    else:
        example = str(e)

    return "# Example: " + example

j2_env.filters["format_descriptions"] = format_descriptions
j2_env.filters["format_comments"] = format_comments
j2_env.filters["required"] = required
j2_env.filters["example"] = example
j2_env.filters["default"] = default

configuration_files = ["config.yaml", "tests/test_data/paired_end/config.yaml",
                       "tests/test_data/dna_repair/config.yaml",
                       "tests/test_data/keep_the_tips/config.yaml",
                       "example/config.yaml"]

base_config = ordered_load(open("utils/config_description.yaml"))
for configuration_file in configuration_files:

    existing_config = yaml.load(open(configuration_file))
    template_config = ordered_load(open("utils/config_description.yaml"))

    for k, v in base_config.items():
        for k2, v2 in list(v.items())[1:]:
            if existing_config.get(k2, ""):
                template_config[k][k2]["default"] = existing_config[k2]

    updated_config_file = j2_env.get_template('utils/base_template.conf').render(
        config=template_config.items()
    )

    time_str = ctime().replace(" ", "_").replace(":", "_")
    subprocess.call("cp {0} {0}_{1}.bkup".format(configuration_file, time_str), shell=True)
    with open(configuration_file, "w+") as outhandle:
        outhandle.write(updated_config_file)
