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


def format_category(desc):
    capitalized = " ".join([s.capitalize() for s in desc.split("_")])
    dashes = "~" * len(desc)
    return " ".join([capitalized, "\n", dashes])

def format_comments(c):

    return "# " + "\n# ".join(wrap(str(c), 78))

def required(r):

    return "# (Required)" if r else "# (Not required)"

def default(d):

    if d != "None" and d is not None:
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

j2_env.filters["format_category"] = format_category
j2_env.filters["format_comments"] = format_comments
j2_env.filters["required"] = required
j2_env.filters["example"] = example
j2_env.filters["default"] = default


base_config = ordered_load(open("utils/config_description.yaml"))
print(base_config.items())

template = j2_env.get_template('docs/configuration_files.j2').render(
    config=base_config.items()
)

open("docs/configuration_files.rst", "w+").write(template)
