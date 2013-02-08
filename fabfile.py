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

import StringIO
import os

from fabric.api import cd, env, execute, put, require, run, task
from fabric.contrib.files import upload_template
from fabric.contrib.project import rsync_project

import diabric.venv
import diabric.config
import diabric.files


HERE = os.path.abspath(os.path.dirname(__file__))


###############
# CONFIGURATION

def post_config(config):
    '''
    Called by one of the deployment environment configuration tasks: dev, prod,
    etc.  Sets some values in `config`.
    '''
    config.venv = os.path.join(config.deploy, 'venv')
    config.data = os.path.join(config.deploy, 'data')
    config.code = os.path.join(config.deploy, 'python')
    config.python = os.path.join(config.venv, 'bin', 'python')
    config.requirements = os.path.join(HERE, 'requirements.txt')
    env.configured = True
    return config


# a global configuration "dictionary" (actually an attribute namespace).
# an alternative to fabric.api.env, which IMHO is a bad place to store
# application configuration, since that can conflict with fabric internals and
# it is not flexible enough to store per-host configuration.
config = diabric.config.Namespace()



@task
def local():
    '''
    Configure tasks for local deployment.
    '''
    env.hosts = ['localhost']
    config.deploy_env = 'local'
    config.deploy = os.path.expanduser('~/deploy/vanvactor_mirna')
    config.system_python = '/usr/local/bin/python'
    post_config(config)


############
# VENV TASKS

@task
def venv_create():
    require('configured')
    diabric.venv.create(config.venv, config.system_python)


@task
def venv_install():
    require('configured')
    diabric.venv.install(config.venv, config.requirements)


@task
def venv_freeze():
    require('configured')
    diabric.venv.freeze(config.venv, config.requirements)


@task
def venv_remove():
    require('configured')
    diabric.venv.remove(config.venv)


############
# DEPLOYMENT

@task
def clean():
    require('configured')
    run('rm -rf ' + config.code)


@task
def init():
    require('configured')
    run('mkdir -p -m 2775 ' + config.data + ' ' + config.code)



@task
def conf():
    '''
    '''
    require('configured')

    # 
    out = StringIO.StringIO()
    out.write("# Autogenerated file.  Do not modify.\n")
    out.write("data_dir = '{}'\n".format(config.data))
    put(out, os.path.join(config.code, 'deployenv.py'), mode=0664)


@task
def deploy():
    '''
    '''
    require('configured')
    put(os.path.join(HERE, 'python'), os.path.join(config.deploy), mode=0664)
    # do not deploy these files/dirs.
    # deployExcludes = ['**/old', '**/old/**', '**/semantic.cache', '**/.svn', '**/.svn/**', '**/*.pyc',
                      # '**/*.pyo', '**/.DS_Store', '**/*~']
    # rsync_project(
        # local_dir=os.path.join(HERE, 'python'),
        # remote_dir=config.deploy,
        # exclude=deployExcludes,
        # delete=False)



@task
def all():
    require('configured')
    execute(clean)
    execute(init)
    execute(conf)
    execute(deploy)


if __name__ == '__main__':
    main()


# pass

