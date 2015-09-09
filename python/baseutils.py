# ----------------------------------------------------------------------------------------
# utilities with minimal dependencies (in principle I should split up utils.py anyway, but a.t.m. I only care that scons and sklearn break each other)

# ----------------------------------------------------------------------------------------
def get_extra_str(extra_list):
    if len(extra_list) == 0:
        return ''
    modified_list = [str(ex).replace(':', ',').replace('--', '__').replace(' ', '+') for ex in extra_list]
    return ' --extra-args ' + ':'.join(modified_list)
