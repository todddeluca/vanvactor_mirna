'''
# Usage examples:

Initialize local deployment environment:

    fab local init

Create virtual environment and install packages:

    fab local venv_create venv_install

Deploy code:

    fab local all

Example of running code:

    cd /Users/td23/deploy/vanvactor_mirna/python
    ../venv/bin/python download_microcosm.py

'''

from __future__ import print_function
import StringIO
import os
import sys

from fabric.api import env, execute, put, require, run, sudo, task, get
import fabric.contrib.files

HERE = os.path.abspath(os.path.dirname(__file__))

import fabvenv
sys.path.append(os.path.join(HERE, 'secrets'))
import prod as secrets



# Stardog database data and license key live in STARDOG_HOME
stardog_version = '1.2.3'
stardog_name = 'stardog-{}'.format(stardog_version)
stardog_home = '/home/ubuntu/data/stardog'


###############
# CONFIGURATION

# 'config' is an alterantive to fabric.api.env, which does not risk having
# application configuration conflict with fabric internals, like 'port'.
class Namespace(object):
    ''' An iterable attribute namespace '''
    def __iter__(self):
        return iter(self.__dict__)

config = Namespace()


@task
def local():
    '''
    Configure tasks for local deployment.
    '''
    env.hosts = ['localhost']
    config.deploy_env = 'local'
    config.proj = os.path.expanduser('~/deploy/vanvactor_mirna')
    config.system_python = '/usr/local/bin/python'
    post_config(config)


@task
def prod():
    '''
    Configure for deployment to EC2 production semantic database server
    '''
    env.hosts = ['ec2-23-21-187-71.compute-1.amazonaws.com']
    env.user = 'ubuntu'
    env.key_filename = os.path.expanduser('~/.ssh/tfd_20120531.pem')
    config.deploy_env = 'prod'
    config.proj = os.path.expanduser('/home/ubuntu/deploy/vanvactor_mirna')
    config.system_python = '/usr/bin/python'
    post_config(config)


def post_config(config):
    '''
    Called by one of the deployment environment configuration tasks: dev, prod,
    etc.  Sets some values in `config`.
    '''
    config.venv = os.path.join(config.proj, 'venv')
    config.data = os.path.join(config.proj, 'data')
    config.code = os.path.join(config.proj, 'python')
    config.python = os.path.join(config.venv, 'bin', 'python')
    config.requirements = os.path.join(HERE, 'requirements.txt')
    config.virtuoso_load_dir = '/home/ubuntu/tmp/virtuoso_load_dir'
    env.configured = True
    return config


###########################
# VIRTUAL ENVIRONMENT TASKS

@task
def venv_create():
    require('configured')
    venv = fabvenv.Venv(config.venv, config.requirements)
    if not venv.exists():
        venv.create(config.system_python)


@task
def venv_install():
    require('configured')
    fabvenv.Venv(config.venv, config.requirements).install()


@task
def venv_upgrade():
    require('configured')
    fabvenv.Venv(config.venv, config.requirements).upgrade()


@task
def venv_freeze():
    require('configured')
    fabvenv.Venv(config.venv, config.requirements).freeze()


@task
def venv_remove():
    require('configured')
    venv = fabvenv.Venv(config.venv, config.requirements)
    if venv.exists():
        venv.remove()


@task
def venv_pth():
    '''
    Add the code directory to the virtualenv sys.path.
    '''
    require('configured')
    fabvenv.Venv(config.venv, config.requirements).venv_pth([config.code])


###################
# HOST SERVER TASKS
# Installing, configuring, starting, and stopping servers on a host

@task
def install_stardog():
    '''
    '''
    # STARDOG_HOME tells stardog where its data files, etc., are.
    fabric.contrib.files.append(
        '~/.bash_profile', 'export STARDOG_HOME="{}"'.format(stardog_home))

    sd = stardog_name

    # Upload and unzip stardog archive
    sudo('rm -rf ~/{sd}'.format(sd=sd))
    put('installs/{sd}/{sd}.zip'.format(sd=sd), '~')
    run('unzip {sd}.zip'.format(sd=sd))

    # Move stardog dir to /usr/local and link executables to /usr/local/bin
    # so they are on PATH
    uninstall_stardog()
    sudo('mv ~/{sd} /usr/local'.format(sd=sd))
    sudo('ln -s /usr/local/{sd}/stardog /usr/local/{sd}/stardog-admin /usr/local/bin'.format(sd=sd))

    # Create STARDOG_HOME dir for stardog data
    run('mkdir -p {}'.format(stardog_home))
    # Add license key to STARDOG_HOME so stardog will work.
    put('installs/{sd}/stardog-license-key.bin'.format(sd=sd), stardog_home)

    # Change default password for default admin user to something not published
    # on the web
    # This starts stardog as a background process, changes the password, then
    # stops stardog (using then new password).  All of this needs to be done
    # in the same run command (except the stop?) b/c stardog starts as a 
    # background process which dies when the run command ends.
    run('''stardog-admin server start --username admin --passwd admin && \\
           stardog-admin user passwd --username admin --passwd admin \\
           --new-password {password} && \\
           stardog-admin server stop --username admin --passwd {password} \\
        '''.format(password=secrets.stardog_admin_password))

@task
def conf_stardog():
    put(os.path.join(HERE, 'conf', stardog_name, 'stardog.properties'),
        stardog_home)

@task
def uninstall_stardog():
    '''
    Remove stardog binaries, leaving stardog data files intact.
    '''
    sudo('rm -rf /usr/local/{sd}'.format(sd=stardog_name))
    sudo('rm -f /usr/local/bin/stardog /usr/local/bin/stardog-admin')

@task
def uninstall_stardog_home():
    '''
    Remove STARDOG_HOME dir, including all stardog database files!
    '''
    sudo('rm -rf {}'.format(stardog_home))



##########################
# RELEASE DEPLOYMENT TASKS


@task
def clean():
    require('configured')
    run('rm -rf ' + config.code)


@task
def init():
    require('configured')
    dirs = [config.data, config.code]
    run('mkdir -p -m 2775 ' + ' '.join(dirs))


@task
def conf():
    '''
    '''
    require('configured')

    # Copy secrets files
    put(os.path.join(HERE, 'secrets/{}.py'.format(config.deploy_env)),
        os.path.join(config.code, 'secrets.py'), mode=0660)

    # Generate deployenv.py
    out = StringIO.StringIO()
    out.write("# Autogenerated file.  Do not modify.\n")
    out.write("datadir = '{}'\n".format(config.data))
    out.write("virtuoso_load_dir = '{}'\n".format(config.virtuoso_load_dir))
    put(out, os.path.join(config.code, 'deployenv.py'), mode=0664)


@task
def deploy():
    '''
    '''
    require('configured')
    put(os.path.join(HERE, 'python'), os.path.join(config.proj), mode=0664)


#########################################
# CONVENIENCE TASKS FOR DOING MANY THINGS


@task()
def copy_results_to_dropbox(results_dir):
    get(os.path.join(config.data, results_dir),
        os.path.expanduser('~/Dropbox/MIRConservation'))


@task
def full():
    require('configured')
    execute(venv_remove)
    execute(venv_create)
    execute(most)


@task
def most():
    require('configured')
    execute(venv_install)
    execute(venv_pth)
    execute(clean)
    execute(init)
    execute(conf)
    execute(deploy)


