from slitlessutils import config

cfg = config.Config()

# download the latest reference files
reffile = cfg.retrieve_reffiles(update=True)
