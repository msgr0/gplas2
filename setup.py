from setuptools import setup, find_packages

with open('envs/requirements.txt') as file_open:
     requirements = file_open.read().splitlines()

setup(
    name="gplas",
    version="1.1.2",
    scripts=["gplas/snakefiles/mlplasmidssnake.smk"],
    packages=find_packages(),
    install_requires=requirements,
    include_package_data=True,
    package_data={'': ['gplas/*']}

)


