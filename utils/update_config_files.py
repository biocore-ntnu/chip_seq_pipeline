
from jinja2 import Environment, FileSystemLoader

env = Environment(loader=FileSystemLoader("utils/"))
env.get_template("config_template.yaml").render()
