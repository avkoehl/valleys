import sys

import toml

from valleys.workflow import full_workflow

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python main.py <config_file>')
        sys.exit(1)

    CONFIG_FILE = sys.argv[1]

    config = toml.load(CONFIG_FILE)

    full_workflow(config)
