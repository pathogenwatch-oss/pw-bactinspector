from subprocess import Popen, PIPE
import os
def run_command(command_and_arguments, shell = False):
    """
    run a command with subprocess having supplied commands in a list named command_and_arguments
    """
    proc = Popen(command_and_arguments,
        stdout = PIPE,
        stderr = PIPE,
        shell = shell
    )
    stdout, stderr = proc.communicate()
 
    return proc.returncode, stdout, stderr

def get_base_name(filename):
    """
    get the basename for a file without extension
    """
    return '.'.join(os.path.basename(filename).split('.')[:-1])

def add_new_file_extension(filename, new_extension):
    """
    remove extension and add another
    """
    basename = get_base_name(filename)
    dirname = os.path.dirname(filename)
    return os.path.join(dirname, '{0}.{1}'.format(basename, new_extension))
