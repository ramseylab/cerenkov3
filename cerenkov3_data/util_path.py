import os

_cur_dir = os.path.dirname(os.path.realpath(__file__))

def get_path(sub_dir):
    return os.path.join(_cur_dir, sub_dir)