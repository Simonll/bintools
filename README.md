# run

## install in conda env
```bash
pip install git+https://github.com/Simonll/bintools.git
```
## pre-commit install
- https://pre-commit.com/
```bash
conda install -c conda-forge pre-commit
```
- copy .pre-commit-config.yaml from another repository
- copy pyproject.toml form another repository
```bash
pre-commit install
```
## create conda environment
- cheat sheet https://docs.conda.io/projects/conda/en/latest/_downloads/843d9e0198f2a193a3484886fa28163c/conda-cheatsheet.pdf
```bash
conda env create --file environment.yml
```
## pruning cache from docker
```bash
docker builder prune
docker system prune -a
```
