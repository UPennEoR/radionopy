from setuptools import find_packages
__version__ = '1.0.1-dev'

setup_args = {
    'name': 'radionopy',
    'author': 'Immanuel Washington',
    'author_email': 'immwa at sas.upenn.edu',
    'description': 'package for ionosphere RM',
    'url': 'https://github.com/jaguirre/radionopy.git',
    'license': '?',
    'package_dir' : {'radionopy': ''},
    'packages' : find_packages(),
    #'install_requires': ['alembic>=0.8.2', 'Flask>=0.10.1', 'Flask-Login>=0.2.11', 'Flask-Migrate>=1.5.1', 'Flask-Script>=2.0.5',
    #                       'Flask-SQLAlchemy>=2.0', 'itsdangerous>=0.24', 'Jinja2>=2.8', 'Mako>=1.0.2', 'MarkupSafe>=0.23',
    #                       'mysqlclient>=1.3.6', 'psycopg2>=2.6.1', 'python-editor>=0.4', 'requests>=2.7.0', 'requests-futures>=0.9.5',
    #                       'SQLAlchemy>=1.0.8', 'Werkzeug>=0.10.4', 'wheel>=0.24.0', 'paramiko', 'prettytable', 'numpy'],
    'version': __version__,
}

if __name__ == '__main__':
    try:
        from setuptools import setup
    except:
        from distutils.core import setup
    try:
        apply(setup, (), setup_args)
    except:
        setup(**setup_args)
