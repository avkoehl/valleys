import toml

from valleys.workflow import full_workflow

CONFIG_FILE = './configs/1801010701.toml'

config = toml.load(CONFIG_FILE)

full_workflow(config)
